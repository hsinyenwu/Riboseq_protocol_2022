You can run your code on locally (PC/Mac), linux server or cluster with properly installed/loaded software packages.
You have to decide the $INPUT and $OUTPUT for each step.

### Prepare example file from the sequencing file
The NEB1.fastq.gz is available on NCBI SRA BioProject ID PRJNA854638 after the project is published.
```
zcat NEB1.fastq.gz | head -1000000 > NEB1.h1M.fastq.gz #Extract first 1 million lines
```

### Step 1 code (FASTQC for quality check):
```
#$OUTPUT output directory
fastqc -o $OUTPUT -t 10 NEB1.h1M.fastq.gz
```

### Step 2 code (FASTX toolkit to trim adaptor sequence):
```
#$INPUT input directory
#$OUTPUT output directory
#Please see FASTX toolkit online page for the parameters used here 
zcat $INPUT/NEB1.h1M.fastq.gz | fastx_clipper -o $OUTPUT/NEB1.h1M.trim.fastq -a CTGTAGGCACCATCAAT -l 20 -c -n -v
```

### Step 3a code (Bowtie2: Build index for contamination sequences)
```
#$CF is the directory contain FASTA files for contamination sequence 
#$OUTPUT output directory
cd $OUTPUT
bowtie2-build $CF/Araport11_contam5.fa Contam5
```
### Step 3b code (Bowtie2: remove contamination sequences)
```
#$Contam5 is the directory with the index for the contamination sequences 
bowtie2 -L 20 -p 8 -x $Contam5 $INPUT/NEB1.h1M.trim.fastq --un-gz $OUTPUT/NEB1.h1M.trim.noContam5.fastq.gz
```

### Step 4a code (STAR: create genome/transcriptome index for STAR)
```
#$starIndex is the directory contain FASTA files
#$FASTA is the path to the TAIR10 genome FASTA file
#$GTF is the path to the Araport11 GTF file 
#$OUTPUT is the output directory
#Please see STAR aligner online page for the parameters used here 

cd $OUTPUT
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $starIndex \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 34 \
```

### Step 4b code (STAR mapping)
```
#$starIndex is the directory contain index files
#$INPUT is the directory with the NEB1.h1M.trim.noContam5.fastq.gz file
#$OUTPUT is the output directory

STAR --runThreadN 10 \
--genomeDir $starIndex \
--readFilesCommand zcat \
--readFilesIn $INPUT/NEB1.h1M.trim.noContam5.fastq.gz \
--alignIntronMax 5000 \
--alignIntronMin 15 \
--outFilterMismatchNmax 1 \
--outFilterMultimapNmax 20 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 2 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--outFileNamePrefix "star_ribo_NEB1_h1M" \
```












