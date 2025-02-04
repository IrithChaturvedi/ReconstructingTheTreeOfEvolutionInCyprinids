# Reconstructing The Tree Of Evolution In Cyprinids

## Identifying the genes that were targeted for sequencing in our previous study on Cyprinidae

### PROJECT 1
#### Identifying the genes that were targeted for sequencing in our previous study:

### Setup Box and VPN:
- **Box:** [Box Login](https://uofi.app.box.com/)
- **VPN:** [Illinois VPN Setup](https://answers.uillinois.edu/illinois/98773)
- **Login to ABE server:** `ssh <netid>@abe.life.illinois.edu` (Use VPN for off-campus access)

### Install Necessary Software (For WSL on ABE)
#### If using Windows:
- **Install Windows Subsystem for Linux (WSL):** [YouTube Tutorial](https://www.youtube.com/watch?v=IgMnBSW_aUM&authuser=0)

#### Software Installations:
- **Miniconda:** [Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
- **goalign** (works with FASTA alignments):
  ```sh
  conda install -c bioconda goalign
  ```
- **seqkit** (works with FASTQ):
  ```sh
  conda install -c bioconda seqkit
  ```
- **fastq-dl** (downloads FASTQ files from SRA database):
  - [GitHub Repository](https://github.com/rpetit3/fastq-dl)
  ```sh
  conda install -c bioconda fastq-dl
  ```
- **Trimmomatic**:
  - [Trimmomatic Website](http://www.usadellab.org/cms/?page=trimmomatic)
  ```sh
  conda install -c bioconda trimmomatic
  ```
- **SPAdes:**
  - [SPAdes Repository](https://github.com/ablab/spades)
  ```sh
  conda install -c bioconda spades
  ```
- **CDhit** (removes duplicate sequences):
  - [CDhit Documentation](https://bioconda.github.io/recipes/cd-hit/README.html)
  ```sh
  conda install -c bioconda cd-hit
  ```
- **NCBI BLAST** (matches sequences by similarity):
  - [BLAST on Anaconda](https://anaconda.org/bioconda/blast)
  ```sh
  conda install -c bioconda blast
  ```

---

## Step-by-Step Process

### STEP 1: Download Danio tinwini dataset
- **Dataset ID:** SRR4556492 ([NCBI Link](https://www.ncbi.nlm.nih.gov/sra/SRX2309841))
- Use `fastq-dl` to download the files:
  ```sh
  conda activate fastq-dl
  fastq-dl --accession <number> --provider sra
  ```

### STEP 2: Trim Danio tinwini data using Trimmomatic
```sh
trimmomatic PE -threads 4 -phred33 -trimlog sample.trimlog \
SRR###_R1.fastq SRR###_R2.fastq \
SRR###_R1.trimmed.fastq SRR###_rem1.fastq.gz \
SRR###_R2.trimmed.fastq.gz SRR###_rem2.fastq.gz \
ILLUMINACLIP:adapters.fa:2:30:10 LEADING:5 TRAILING:5 \
SLIDINGWINDOW:4:15 MINLEN:31
```

### STEP 3: Assemble Danio tinwini data using SPAdes
```sh
source activate py27
nohup spades.py --pe1-1 <trimmed fastq file(1)> --pe1-2 <trimmed fastq file(2)> -o <desired output file> &
```
- **Output Files:** `scaffolds.fasta`, `contigs.fasta`

### STEP 4: Use CD-HIT to remove redundant sequences
```sh
cd-hit -i input.fasta -o output.cdhit -c 1
```

### STEP 5: BLAST Danio tinwini assembled sequences to Danio rerio genome
- **Download the Danio rerio genome:** [ENSEMBL](https://useast.ensembl.org/Danio_rerio/Info/Index)
- **Run BLAST:**
```sh
makeblastdb -dbtype nucl -in Danio_rerio_genome.fasta -out Danio_rerio_genome
blastn -db Danio_rerio_genome -query Danio_tinwini_asm_dedup.fasta -outfmt 6 -out Danio_tinwini_asm_dedup_V_Danio_rerio_genome_blastn.txt
```

### STEP 6: Compare BLAST results with Danio rerio GFF3 annotations
```sh
awk '{print($2"\t"$9"\t"$10"\t"$1)}' <BLAST output> > <final output file>
bedtools intersect -a blast_results.bed -b annotation.gff3 -wa -wb > intersected_hits.bed
```

### STEP 7: Tab-delimit last column for ENSDARG & ENSDARE analysis
- **Script Provided:** `project_1.py`
```sh
project_1.intersected_hits_tab_delimited(lines)
```

### STEP 8: Repeat Steps 1-7 for other Cyprinid species
- **Species Used:**
  - Danio Feegradei - SRR4556469
  - Danio Margaritatus - SRR4556491
  - Puntius Sophore - SRR4578774

### STEP 9: Find intersect with Danio Lemmon sequences
```sh
nohup blastn -db Danio_rerio_primaryassembly.db -query Danio_Lemmon_sequences -outfmt 6 -out Danio_Lemmon_asm_dedup_V_Danio_rerio_genome_blastn.txt &
```

### STEP 10: Intersect hits from different species with Danio Lemmon sequences
```sh
bedtools intersect -a intersected_hits_DANIO.bed -b intersected_hits_Lemmon.bed -wa -wb > intersected_hits_DANIO_and_Lemmon.bed
```

### STEP 11: Analyze resultant files
- **Plot Coverage vs Length:**
  ```sh
  project_1.make_coverage_dataframe_genes(lines, feature_index, ensembl_id_index)
  project_1.make_coverage_dataframe_exons(lines, feature_index, ensembl_id_index)
  ```

### STEP 12: Parse files for 'good' coverage range
```sh
project_1.good_coverage_range(lines, high_cov, low_cov)
```

### STEP 13: Remove redundant gene and exon features
```sh
project_1.redundant_gene_exon_parse(lines, feature_sort, feature_index)
```

### STEP 14: Additional Analysis
- Create a **binary presence-absence table** to analyze unique features across species.

---

### Summary
This project reconstructs the **Tree of Evolution in Cyprinids** by performing de-novo assembly and sequence alignment of various species. Using bioinformatics tools like **SPAdes, CD-HIT, and BLAST**, the study identifies genes and exon features conserved across species, enabling further evolutionary analysis.

**Scripts Provided:** `project_1.py`, `BINARY_TABLE script`
