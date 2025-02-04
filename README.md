# Reconstructing The Tree Of Evolution In Cyprinids
Identifying the genes that were targeted for sequencing in our previous study on Cyprinidae
PROJECT 1
Identifying the genes that were targeted for sequencing
in our previous study:

Setup Box and VPN:

Box - https://uofi.app.box.com/
VPN - https://answers.uillinois.edu/illinois/98773
Login to ABE server, ssh <netid>@abe.life.illinois.edu, Use VPN for off campus access

Install all necessary software - For WSL on abe

If using Windows - Install Windows Subsystem for Linux (WSL) - https://www.youtube.com/watch?v=IgMnBSW_aUM&authuser=0

Miniconda - https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

goalign - works with FASTA alignments, needs a flag for unaligned sequences

conda install -c bioconda goalign

seqkit - works with FASTQ

conda install -c bioconda seqkit

fastq-dl - allow download FASTQ files from SRA database: https://github.com/rpetit3/fastq-dl

conda install -c bioconda fastq-dl

Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
Use adapters.fa with trimmomatic

conda install -c bioconda trimmomatic
SPAdes: https://github.com/ablab/spades

conda install -c bioconda spades

CDhit: remove duplicate sequences
https://bioconda.github.io/recipes/cd-hit/README.html

conda install -c bioconda cd-hit

ncbi BLAST: Matching sequences by similarity:
https://anaconda.org/bioconda/blast

conda install -c bioconda blast

STEP - 1
Download Danio tinwini dataset: SRR4556492
https://www.ncbi.nlm.nih.gov/sra/SRX2309841[accn]
Use fastq-dl to download the files
Goal will be two fastq.gz files, one for each end
conda activate fastq-dl
fastq-dl --accession <number> --provider sra

STEP - 2
Trim Danio tinwini data using Trimmomatic
if installed using conda rather than downloading the jar file, then the command should be trimmomatic instead
Given a sample with name “sample” and fastq files starting with some SRR number:

trimmomatic PE -threads 4 -phred33 -trimlog sample.trimlog SRR###_R1.fastq SRR###_R2.fastq SRR###_R1.trimmed.fastq SRR###_rem1.fastq.gz SRR###_R2.trimmed.fastq.gz SRR###_rem2.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:31

Expected - 4 files 2 trimmed and 2 rem files

STEP - 3 (Use python 2.7 as base to install SPAdes)
Assemble Danio tinwini data using SPAdes
If possible, keep track of runtime
Expected result should be a FASTA file containing numerous (thousands) of sequence headers (beginning with >). You may also want to store a log file.
SPAdes takes a while, so you might want to consider using something like screen, byobu, or nohup to run this in the background

source activate py27

nohup <command> &- Lets long processes run without interruptions in the background

nohup spades.py --pe1-1 <trimmed fastq file(1)> --pe1-2 <trimmed fastq file(2)> -o <desired output file> &

scaffolds.fasta and contigs.fasta are produced

STEP - 4 (nohup to be used)
Use CD-HIT to identify and remove redundant sequences from assembled Danio tinwini data

cd-hit -i input.fasta -o output.cdhit -c 1

cdhit should be a fasta file, technically

STEP - 5
BLAST Danio tinwini assembled, non-redundant sequences to Danio rerio genome
Download the Danio rerio genome: https://useast.ensembl.org/Danio_rerio/Info/Index
Download the dna_rm.primary_assembly genome (fasta) as well as the annotation (GFF3)
annotation: Danio_rerio.GRCz11.110.gff3.gz
Tutorial: https://conmeehan.github.io/blast+tutorial.html
makeblastdb on the Danio rerio genome (FASTA) file

makeblastdb -dbtype nucl -in Danio_rerio_genome.fasta -out Danio_rerio_genome

blastn -db Danio_rerio_genome -query Danio_tinwini_asm_dedup.fasta -outfmt 6 -out Danio_tinwini_asm_dedup_V_Danio_rerio_genome_blastn.txt
nohup blastn -db Danio_rerio_primaryassembly.db -query cd-hit_danio_res -outfmt 6 -out Danio_tinwini_asm_dedup_V_Danio_rerio_genome_blastn.txt &

-db: replace with danio rerio genome blast database

Use blastn with D tinwini deduplicated assembly as query and Danio rerio genome database
Use -outfmt 6

STEP - 6
Compare output table to Danio rerio GFF3 annotations to identify contigs

Using BLAST tabular results, compare to Danio rerio GFF3 using bedtools
Convert BLAST result to bed format: GitHub - nterhoeven/blast2bed: convert a blast output to a bed file

Parse BLAST output file - 2nd and 3rd columns should be integers (length)
For bedtools to work, file should be TAB delimited

command: awk ‘{print($2"\t"$9"\t"$10"\t"$1)}’ <BLAST output> > <final output file>

Sort the 2nd and 3rd column - second column data to be less then third column for the coordinates
(script provided)

Also parse BLAST output to contain only those alignments with e-values less than 1e-5. (script provided)

Parse .gff3 database to contain only gene and exon features. (script ReadMe provided), Function:
project_1.gff3_sort_key(lines)

bedtools intersect -a blast_results.bed -b annotation.gff3 -wa -wb > intersected_hits.bed

STEP - 7
Last column needs to be tab delimited to analyze ENSDARG and ENSDARE numbers for exons and genes. Script provided. Check module project_1.py following function:
project_1.intersected_hits_tab_delimited(lines)

STEP - 8
Follow Steps 1 to 7 for different species to conduct De-novo assembly of genes and exon features targeted in a previous assembly.
In this case, species of the same family, here Cyprinidae are chosen, with the hypothesis that similar enough species will usually contain the same gene features. See functions in project_1.py module.

Species Used:
    Danio Feegradei - SRR4556469 
    Danio Margaritatus - SRR4556491
    Puntius Sophore - SRR4578774
    (download all following STEP 1)
STEP - 9
Using a reference set of genes which were already found to be good enough to be targeted. Here the Danio Lemmon data, created by Alan and Emily Lemmon for the Danio rerio species is used to find an intersect between the species to find common features.
Creating the Danio Lemmon intersected hits file.
The Lemmon sequences are presented as a FASTA file. The sequences are to be BLAST against the Danio rerio gff3 file to create an intersected hits file for the same.
Follow STEP 5 to 7, BLAST command example:
nohup blastn -db Danio_rerio_primaryassembly.db -query Danio_Lemmon_sequences -outfmt 6 -out Danio_Lemmon_asm_dedup_V_Danio_rerio_genome_blastn.txt &

STEP - 10
The intersected hits file of all species are to be INTERSECTED AGAIN with the Danio Lemmon intersected hits file, to find a good list of sequences with features that may have been targeted. Example:
bedtools intersect -a intersected_hits_DANIO.bed -b intersected_hits_Lemmon.bed -wa -wb > intersected_hits_DANIO_and_Lemmon.bed

STEP - 11
The resultant files are then analyzed.
Conversion of columns delimited with ‘_’ (underscore) to columns delimited with space. Script provided, Check module project_1.py following function:

Plotting coverage vs length graph to figure out the coverage range to be selected for the features of each species. The region where the graph is biased or highly concentrated is to be selected as the coverage range for that species. For easier analysis, it may also make sense to create a dataframe for each file containing only the most important features.
    Script provided, Check module project_1.py following function:
    project_1.make_coverage_dataframe_genes(lines, feature_index, ensembl_id_index) - GENES

    project_1.make_coverage_dataframe_exons(lines, feature_index, ensembl_id_index) - EXONS
STEP - 12
Parsing files to only contain sequences in the ‘good’ coverage range resultant from STEP 11.
Script provided, check module project_1.py function - project_1.good_coverage_range(lines, high_cov, low_cov)

STEP - 13
Removal of redundant gene and exon features.
Redundant gene and exon features present within the same coordinates is necessary to produce a much finished and produced result. Simply pass the files through the module function provided in the script.
Function: (Check script read_me for specifications)
project_1.redundant_gene_exon_parse(lines, feature_sort, feature_index)

STEP - 14
Additional Analysis:
To perform additional analysis and find the unique features present in two or more species, A binary presence-absence table can be created by passing the good intersected hits files through the following methods provided in the BINARY TABLE script.
