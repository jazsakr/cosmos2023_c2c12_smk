import pandas as pd

# read sample file and get sample names
metadata = pd.read_csv('c2c12_samples.csv', sep=',')
metadata_df = pd.DataFrame(data = metadata)
sn_list = metadata_df['sample_file_name'].tolist()
exp_list = list()
for line in sn_list:
	exp_list.append(line.split('_')[0])
exp = metadata_df['experiment'].tolist()

# set reference names:
ref_fa_name = config['reference']['mouse']['fa_name'],
ref_gtf_name = config['reference']['mouse']['gtf_name']

# set experiment and talon database name:
talon_db_name = '_talon'
talon_list_name = '_talon-list.csv'

rule all:
	input:
		expand(config['talon']['talon_ab_f'], sample=sn_list, exp=exp)

# map the reads to reference genome using minimap2 and parameters for nanopore reads

rule minimap:
	input: 
		fastq = config['sample']['fastq'],
		ref_fa = config['reference']['mouse']['fa'],
		ref_bed = config['reference']['mouse']['bed']
	output: 
		sam = config['sample']['sam']
	threads: 16
	resources:
		mem_gb = 32
	log: 'data/mapped/logs/minimap_{sample}.log'
	shell:
		'minimap2 --MD -t {threads} -ax splice --junc-bed {input.ref_bed} \
		{input.ref_fa} {input.fastq} > {output.sam} 2> {log}'

# convert sam to bam file using samtools

rule sam_to_bam:
	input: 
		sam = config['sample']['sam']
	output:
		bam = config['sample']['bam']
	threads: 4
	resources:
		mem_gb = 16
	run: 
		shell('samtools sort --threads={threads} -O bam {input.sam} > {output.bam} \
		&& samtools index -@ {threads} {output.bam}')

# reverse flipped reads with script

rule sam_reverse:
	input:
		bam = config['sample']['bam']
	output:
		sam_reverse = config['sample']['sam_reverse']
	threads: 4
	resources:
		mem_gb = 16
	script: 
		'scripts/reverse_bam_smk.py'

# reference-based error correction by removing microindels, mismatches, and 
# noncanonical splice junctions in a variant-aware manner using TranscriptClean

rule transcriptclean:
	input:
		sam_reverse = config['sample']['sam_reverse'],
		tc = config['packages']['transcriptclean'],
		ref_fa = config['reference']['mouse']['fa'],
		ref_sj = config['reference']['mouse']['sj']
	output:
		clean = config['sample']['clean']
	params:
		tc_opref = config['sample']['clean'].rsplit('_clean.sam', maxsplit=1)[0]
	threads: 16
	resources:
		mem_gb = 80
	shell: 
		'mkdir {params.tc_opref}_tc-tmp && \
		python {input.tc}/TranscriptClean.py \
		-t {threads} \
		--sam {input.sam_reverse} \
		--genome {input.ref_fa} \
		--spliceJns {input.ref_sj} \
		--canonOnly \
		--primaryOnly \
		--deleteTmp \
		--outprefix ./{params.tc_opref}_tc-tmp/{wildcards.sample}r && \
		mv ./{params.tc_opref}_tc-tmp/{wildcards.sample}r* ./data/processed/ \
		&& rm -r ./{params.tc_opref}_tc-tmp'


"""
TALON
"""

# flag reads for internal priming by recording fraction As in a custom fA:f SAM tag

rule talon_label_reads:
	input:
		talon = config['packages']['talon'],
		clean = config['sample']['clean'],
		ref_fa = config['reference']['mouse']['fa']
	output:
		labeled = config['sample']['labeled']
	params:
		talon_opref = config['sample']['labeled'].rsplit('_labeled.sam', maxsplit=1)[0]
	threads: 16
	resources:
		mem_gb = 32
	shell:
		'talon_label_reads \
		--f {input.clean} \
		--g {input.ref_fa} \
		--t {threads} \
		--tmpDir ./{params.talon_opref}_talon-tmp \
		--ar 20  \
		--deleteTmp \
		--o ./{params.talon_opref}'

# create a comma-delimited configuration file for TALON from the sample sheet

rule talon_config:
	input: 
		get_all = expand(config['sample']['labeled'], sample=sn_list)
	output:
		talon_config = config['talon']['talon_config']
	threads: 1
	resources:
		mem_gb = 4
	script:
		'scripts/talon_config_smk.py'

# initialize talon database from GTF annotation

rule talon_db:
	input:
		talon = config['packages']['talon'],
		ref_gtf = config['reference']['mouse']['gtf']
	output:
		db = config['talon']['talon_db']
	params:
		talon_opref = config['talon']['talon_ab'].rsplit('_talon', maxsplit=1)[0]
	threads: 1
	resources:
		mem_gb = 4
	shell: 
		'talon_initialize_database \
		--f {input.ref_gtf} \
		--g {ref_fa_name} \
		--a {ref_gtf_name} \
		--l 0 \
		--idprefix TALONT \
		--5p 500 \
		--3p 300 \
		--o ./{params.talon_opref}{talon_db_name}'

# annotate your reads and assign transcript model in the TALON database

rule talon:
	input:
		talon = config['packages']['talon'],
		talon_db = config['talon']['talon_db'],
		talon_config = config['talon']['talon_config']
	output:
		talon_ra = config['talon']['talon_read_annot']
	params:
		talon_opref = config['talon']['talon_ab'].rsplit('_talon', maxsplit=1)[0]
	threads: 16
	resources:
		mem_gb = 120
	shell: 
		'talon \
		--f {input.talon_config} \
		--db {input.talon_db} \
		--build {ref_fa_name} \
		--tmpDir ./{params.talon_opref}_talon-tmp \
		--threads {threads} \
		--o ./{params.talon_opref} \
		&& rm -r ./{params.talon_opref}_talon-tmp'

# extract raw transcript count matrix from TALON database

rule talon_abundance:
	input:
		talon = config['packages']['talon'],
		talon_db = config['talon']['talon_db'],
		talon_ra = config['talon']['talon_read_annot']
	output: 
		talon_ab = config['talon']['talon_ab']
	params:
		talon_opref = config['talon']['talon_ab'].rsplit('_talon', maxsplit=1)[0]
	threads: 1
	resources:
		mem_gb = 16
	shell: 
		'talon_abundance \
		--db {input.talon_db} \
		-a {ref_gtf_name} \
		-b {ref_fa_name} \
		--o ./{params.talon_opref}'

# generates a whitelist to use for filtering the transcript count matrix

rule talon_filter_list:
	input:
		talon_ab = config['talon']['talon_ab'],
		talon = config['packages']['talon'],
		talon_db = config['talon']['talon_db']
	output: 
		talon_list = config['talon']['talon_list']
	params:
		talon_opref = config['talon']['talon_ab'].rsplit('_talon', maxsplit=1)[0]
	threads: 1
	resources:
		mem_gb = 16
	shell:
		'talon_filter_transcripts \
		--db {input.talon_db} \
		-a {ref_gtf_name} \
		--maxFracA=0.5 \
		--minCount=3 \
		--minDatasets=2 \
		--o ./{params.talon_opref}{talon_list_name}'

# extract filtered transcript count matrix from TALON database using the whitelist

rule talon_filtered_abundance:
	input:
		talon = config['packages']['talon'],
		talon_db = config['talon']['talon_db'],
		talon_list = config['talon']['talon_list']
	output: 
		talon_af = config['talon']['talon_ab_f']
	params:
		talon_opref = config['talon']['talon_ab'].rsplit('_talon', maxsplit=1)[0]
	threads: 1
	resources:
		mem_gb = 16
	shell:
		'talon_abundance \
		--db {input.talon_db} \
		-a {ref_gtf_name} \
		-b {ref_fa_name} \
		--whitelist {input.talon_list} \
		--o ./{params.talon_opref}'
