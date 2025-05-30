# Elixir Course, Bari, Italy, 2025
## Elixir Course, Bari, Italy, 2025
### Elixir Course, Bari, Italy, 2025
#### Elixir Course, Bari, Italy, 2025
##### Elixir Course, Bari, Italy, 2025
###### [ !! nome del liunk per la pipeline di Davide !! ](https://github.com/metashot/mag-illumina)

## Hands-on n.1 - Taxonomic and functional profiling using shotgun data
### Topic n.1: Preprocessing
### Step n.1: get into the right place
```
cd /home/user<YOUR USER NAME>
```
### Step n.2: set up anaconda and check whether your environment is visible
```
DON'T INSTALL IT...
##wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
##bash Anaconda3-2024.10-1-Linux-x86_64.sh

path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
source ${path}/activate

conda info --envs
```
### Step n.3: raw data pre-processing on fastq example files "seq_1.fastq.gz" and "seq_2.fastq.gz" from https://github.com/biobakery/biobakery/wiki/kneaddata
```
##conda create -n <trimmomatic> -c bioconda trimmomatic ## DON'T DO IT. WE DID ALREADY
##conda create -n <bowtie2> -c bioconda bowtie2 ## DON'T DO IT. WE DID ALREADY
##conda create -n <samtools> -c bioconda samtools ## DON'T DO IT. WE DID ALREADY

mkdir 1_pre-processing
cd 1_pre-processing

wget https://github.com/biobakery/kneaddata/files/4703820/input.zip
unzip input.zip
cd input
```

### Step n.4: Define variable "s" with the sampleID and run TRIMMOMATIC
```
s="seq"

source ${path}/activate trimmomatic

trimmomatic PE -threads 8 -phred33 -trimlog ${s}_trimmomatic.log ${s}1.fastq ${s}2.fastq \
${s}_filtered_1.fastq ${s}_unpaired_1.fastq ${s}_filtered_2.fastq ${s}_unpaired_2.fastq \
ILLUMINACLIP:${path}/../envs/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75

for i in *.fastq; do echo -ne "${i}\t"; cat "$i" | wc -l; done
```

### Step n. 5: Generate bowtie2 index of the human genome GCF_009914755.1_T2T-CHM13v2.0.fna (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0.fna)
```
human_gen_path="/home/ubuntu/shotgun_course/human_genome/"
conda deactivate
source ${path}/activate bowtie2

##VERSION 4 HOURS LONG:
## mkdir -p ../human_genome/
##bowtie2-build ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0.fna ../human_genome/GCF_009914755.1_T2T-CHM13v2.0 ### DON'T RUN IT! IT TAKES A FEW HOURS TO BE EXECUTED

##VERSION 10 SECONDS LONG
bowtie2 -x ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0 -1 ${s}_filtered_1.fastq -2 ${s}_filtered_2.fastq \
    -S ${s}.sam --very-sensitive-local -p 8

conda deactivate
source ${path}/activate samtools

samtools view -bS ${s}.sam > ${s}.bam
samtools view -b -f 12 -F 256 ${s}.bam > ${s}.bothunmapped.bam
samtools sort -n -m 5G -@ 2 ${s}.bothunmapped.bam -o ${s}.bothunmapped.sorted.bam
samtools fastq ${s}.bothunmapped.sorted.bam -1 >(gzip > ${s}_filtered.final_1.fastq.gz) -2 >(gzip > ${s}_filtered.final_2.fastq.gz) -0 /dev/null -s /dev/null -n
#rm ${s}.sam; rm ${s}.bam; rm ${s}.bothunmapped.bam; rm ${s}.bothunmapped.sorted.bam ### IF YOU WANT TO REMOVE THE INTERMEDIATE FILES

for i in *.gz; do echo -ne "${i}\t"; zcat "$i" | wc -l; done
```
### Did the preprocessing produce the same exact number of reads in R1 and R2 ?

### Topic n.2: MetaPhlAn 4: taxonomic profiling using marker genes
### Step n.1: Setup correct variables, activate environment and navigate to the right folders
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"

conda deactivate
source ${path}/activate

##conda create -n <mpa> -c conda-forge -c bioconda python=3.7 metaphlan=4.1.0 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate mpa

mkdir 2_metaphlan
cd 2_metaphlan
```

### Step n.2: Download metagenomic samples
```
mpa_db="/home/ubuntu/shotgun_course/metaphlan_databases/"
db_version="mpa_vJun23_CHOCOPhlAnSGB_202403"

wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014476-Supragingival_plaque.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014494-Posterior_fornix.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014459-Stool.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014464-Anterior_nares.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014470-Tongue_dorsum.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014472-Buccal_mucosa.fasta.gz

s="SRS014476-Supragingival_plaque"
```

### Step n.3: Let's have a look at the MetaPhlAn parameters
```
metaphlan -h
```

### Step n.4: Run MetaPhlAn 4
```
metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --stat_q 0.1 --nproc 8 --bowtie2db ${mpa_db} --index ${db_version}

s="SRS014494-Posterior_fornix"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --stat_q 0.1 --nproc 8 --bowtie2db ${mpa_db} --index ${db_version}
s="SRS014459-Stool"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --stat_q 0.1 --nproc 8 --bowtie2db ${mpa_db} --index ${db_version}
s="SRS014464-Anterior_nares"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --stat_q 0.1 --nproc 8 --bowtie2db ${mpa_db} --index ${db_version}
s="SRS014470-Tongue_dorsum"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --stat_q 0.1 --nproc 8 --bowtie2db ${mpa_db} --index ${db_version}
s="SRS014472-Buccal_mucosa"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --stat_q 0.1 --nproc 8 --bowtie2db ${mpa_db} --index ${db_version}

merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
```

### Topic n.3: Kraken + Bracken: taxonomic profiling using k-mers

### Step n.1: Check everything is set up and move to the right directory
```
conda deactivate
source ${path}/activate

source ${path}/activate kraken
```

### Step n.2: Create kraken database
```
kraken_db="/home/ubuntu/shotgun_course/kraken2_database_gtdb_r220"

kraken2-build --db kraken2_gtdb_r220 --add-to-library gtdb_sequences.fna
kraken2-build --build --db kraken2_gtdb_r220
```

### Step n.3: Let's have a look at Kraken parameters
```
kraken -h
```

### Step n.4: Kraken database build
```

```

### Step n.5: Run Kraken + Braken


### Topic n.4: HUMAnN 4: functional profiling at the community level
### Step n.1: Get into the right directory
```
conda deactivate
source ${path}/activate

## conda create -n <humann> -c bioconda python=3.9 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate humann
## conda install -c biobakery <humann> ## DON'T DO IT. WE DID ALREADY

mkdir 6_humann
cd 6_humann
```

### Step n.2: test that HUMAnN runs properly and have a look at the HUMAnN parameters
```
humann_test
humann_config

humann -h
```

### Step n.3: get a sample from EBI
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/096/SRR15408396/SRR15408396.fastq.gz
```

### Step n.4: RUN humann
```
s="SRR15408396"

## NOW YOU CAN RUN:
## humann --input ${s}.fastq.gz --output ${s} --threads 8

## BUT IT TAKES THREE HOURS... OR YOU CAN RUN:
mkdir -p ${s}

cp /home/ubuntu/course_backup/course/6_humann/${s}/${s}_genefamilies.tsv ${s}/
cp /home/ubuntu/course_backup/course/6_humann/${s}/${s}_pathabundance.tsv ${s}/
cp /home/ubuntu/course_backup/course/6_humann/${s}/${s}_pathcoverage.tsv ${s}/
```

### Step n.5: Manipulate and normalize HUMAnN output
```
humann_renorm_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_genefamilies-relab.tsv -u relab
humann_renorm_table -i ${s}/${s}_pathabundance.tsv -o ${s}/${s}_pathabundance-relab.tsv -u relab
```

### Step n.6: Regrouping genes to other functional categories
```
humann_regroup_table -i ${s}/${s}_genefamilies-relab.tsv -o ${s}/${s}_rxn-relab.tsv --groups uniref90_rxn
```

### Step n.7: Run HUMAnN on a second sample
```
s="SRR15408398"

## SAME:
## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/098/SRR15408398/SRR15408398.fastq.gz
## humann --input ${s}.fastq.gz --output ${s} --threads 8

## FOR NOW, RUN:
mkdir ${s}

cp /home/ubuntu/course_backup/course/6_humann/${s}/${s}_genefamilies.tsv ${s}/
cp /home/ubuntu/course_backup/course/6_humann/${s}/${s}_pathabundance.tsv ${s}/
cp /home/ubuntu/course_backup/course/6_humann/${s}/${s}_pathcoverage.tsv ${s}/

humann_renorm_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_genefamilies-relab.tsv -u relab
humann_renorm_table -i ${s}/${s}_pathabundance.tsv -o ${s}/${s}_pathabundance-relab.tsv -u relab
```

### Step n.8: Merge together community profiles under different ontologies
```
mkdir -p merged

cp SRR15408396/SRR15408396_genefamilies-relab.tsv merged/
cp SRR15408398/SRR15408398_genefamilies-relab.tsv merged/
cp SRR15408396/SRR15408396_pathabundance-relab.tsv merged/
cp SRR15408398/SRR15408398_pathabundance-relab.tsv merged/
cp SRR15408396/SRR15408396_pathcoverage.tsv merged/
cp SRR15408398/SRR15408398_pathcoverage.tsv merged/

humann_join_tables -i merged -o merged/merged_genefamilies-relab.tsv --file_name genefamilies-relab
humann_join_tables -i merged -o merged/merged_pathabundance-relab.tsv --file_name pathabundance-relab
humann_join_tables -i merged -o merged/merged_pathcoverage.tsv --file_name pathcoverage
```

## Hands-on n.2 - Taxonomic profiling beyond the level of species
### Step n.1:

### Step n.2: Getting example files (6 fastq files) from https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS013951.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS014613.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS019161.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS022137.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS055982.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS064276.fastq.bz2
```

### Step n.3: Running MetaPhlAn 4
### Approach n. 1 ==> Running MetaPhlAn 4 to obtain the .sam files of the marker genes' alignments
```
## mpa_db="/home/ubuntu/shotgun_course/metaphlan_databases/"
## db_version="mpa_vJun23_CHOCOPhlAnSGB_202403"

## s="SRS013951"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --bowtie2db ${mpa_db} --index ${db_version}
## s="SRS014613"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --bowtie2db ${mpa_db} --index ${db_version}
## s="SRS019161"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --bowtie2db ${mpa_db} --index ${db_version}
## s="SRS022137"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --bowtie2db ${mpa_db} --index ${db_version}
## s="SRS055982"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --bowtie2db ${mpa_db} --index ${db_version}
## s="SRS064276"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --bowtie2db ${mpa_db} --index ${db_version}
```

### Approach n. 2 ==> Copy the .sam alignments from a pre-existing repository
```
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS013951.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS014613.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS019161.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS022137.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS055982.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS064276.sam.bz2 .
```

### Step n.4: Extract for each sample the alignments over its markers
```
mpa_database="/home/ubuntu/shotgun_course/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl"
sample2markers.py -i *.sam.bz2 -o ./ -n 8 -d ${mpa_database}
```

### Step n.5: Extract marker genes for a species of interest
```
mkdir -p db_markers
```
### Approach n. 1 ==> run the dedicate command
```
extract_markers.py -c t__SGB1877 -o db_markers/ -d ${mpa_database} ## TOO LONG,
```

### Approach n. 2 ==> Copy the pre-built marker files
```
cp /home/ubuntu/course_backup/course/4_strainphlan/db_markers/t__SGB1877.fna db_markers/
```

### Step n.6: Also include a reference genome ("GCF000273725")
```
mkdir -p reference_genomes
wget -P reference_genomes/ http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/reference_genomes/G000273725.fna.bz2
```

### Step n.7: Let's look the StrainPhlAn params
```
strainphlan -h
```

### Step n.8: Run StrainPhlAn 4
```
mkdir -p strainphlan_output
strainphlan -s *.json.bz2 -m db_markers/t__SGB1877.fna -r reference_genomes/G000273725.fna.bz2 -o strainphlan_output -c t__SGB1877 -n 8 -d ${mpa_database}
```

### Step n.9: Let's visualize it ! 
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/metadata.txt
add_metadata_tree.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre -f metadata.txt -m subjectID --string_to_remove .fastq.bz2

conda deactivate
source ${path}/activate graphlan
${path}/../envs/mpa/bin/plot_tree_graphlan.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre.metadata -m subjectID
```

## Hands-on n.3 - Metagenome assembly and binning
### Approach n. 1: follow the protocol
### Step n.1: check everything is set up, download a sample, and run Megahit
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"

conda deactivate
source ${path}/activate

## conda create -n <megahit> -c bioconda megahit ## DON'T DO IT. WE DID ALREADY
source ${path}/activate megahit

mkdir 7_assembly
cd 7_assembly

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/SRR341725/SRR341725_[12].fastq.gz

megahit -h

s="SRR341725"
## MEGAHIT WILL TAKE A FEW HOURS:
## megahit -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -o ${s}.megahit_asm -t 8

## FOR NOW WE CAN COPY THE RESULTS FROM MEGAHIT
mkdir -p ${s}.megahit_asm/

cp /home/ubuntu/course_backup/course/7_assembly/${s}.megahit_asm/final.contigs.fa  ${s}.megahit_asm/
cp /home/ubuntu/course_backup/course/7_assembly/${s}.megahit_asm/contigs.fasta  ${s}.megahit_asm/

## WE ALSO NEED TWO CUSTOM SCRIPT:
cp /home/ubuntu/course_backup/course/7_assembly/filter_contigs.py .
cp /home/ubuntu/course_backup/course/7_assembly/megahit2spades.py .
```

### Step n.2: Binning, i.e. grouping assemblies into genomes using MetaBat2
```
source ${path}/activate
 
## conda create -n <metabat2> -c bioconda metabat2 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate metabat2

## conda install -c bioconda <bowtie2> ## DON'T DO IT. WE DID ALREADY
## conda install -c bioconda <samtools> ## DON'T DO IT. WE DID ALREADY

mkdir 8_MAG-reconstruction
cd 8_MAG-reconstruction

s="SRR341725"

cp ../7_assembly/SRR341725.megahit_asm/contigs_filtered.fasta ./
cp ../7_assembly/SRR341725_1.fastq.gz ./
cp ../7_assembly/SRR341725_2.fastq.gz ./

bowtie2-build contigs_filtered.fasta contigs_filtered
bowtie2 -x contigs_filtered -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -S ${s}.sam -p 8 2> ${s}.bowtie2.log

### THE FOLLOWING MIGHT RAISE DISK ISSUE:
## samtools view -bS ${s}.sam > ${s}.bam
## samtools sort ${s}.bam -o sorted_${s}.bam

## COPY THE RESULT FOR NOW:
cp /home/ubuntu/course_backup/course/8_MAG-reconstruction/sorted_SRR341725.bam .

jgi_summarize_bam_contig_depths --outputDepth ${s}_depth.txt sorted_${s}.bam 2> ${s}_depth.log
```

### Step n.3: Run MetaBat 2 for binning
```
metabat2 -i contigs_filtered.fasta -a ${s}_depth.txt -o ${s}_bins/bin -m 1500 --unbinned -t 8 > ${s}_metabat2.log
```

### Step n.4: Estimate MAG quality using checkM2
```
conda deactivate

## conda create -n <checkm2> -c bioconda checkm2 ## DON'T DO IT. WE DID ALREADY

source ${path}/activate checkm2
## pip install absl-py==1.1.0 ## DON'T DO IT. WE DID ALREADY

## LET'S NOT DOWNLOAD THE DATABASE
## checkm2 database --download --path ./

## WE CAN USE A COPY
checkm2_db="/home/ubuntu/course_backup/course/8_MAG-reconstruction/CheckM2_database/uniref100.KO.1.dmnd"
checkm2 testrun --database_path ${checkm2_db} --threads 8

checkm2 predict -i SRR341725_bins -o SRR341725_checkm2 -x .fa --database_path ${checkm2_db} --threads 8

awk -F'\t' '$2 > 50 && $3 < 5' SRR341725_checkm2/quality_report.tsv > SRR341725_checkm2/quality_report_filtered.tsv

mkdir -p ${s}_bins_filtered
cut -f1 SRR341725_checkm2/quality_report_filtered.tsv | while read -r value; do cp ${s}_bins/${value}.fa ${s}_bins_filtered/; done
```

## Approach n. 2: run the nextflow workflow

Warning: the nextflow commands are here for explanatory purposes only: each of them will require from hours to days to run and considerable computational resouces. The results have been pre-computed and are already present on your workspace.

Make sure that you have nextflow and singularity installed

https://www.nextflow.io/

https://docs.sylabs.io/guides/latest/user-guide/


### Step 1: Assembly and binning


```
cd mag
 nextflow run metashot/mag-illumina \
   --assembler "megahit" \
   --reads '../fastq/*_R{1,2}.fastq.gz' \
   --outdir results \ 
   -with-report report.html
```
This commamnd runs the assembly (using megahit) and binning workflow using pired-ends read files in the ../fastq/ directory. 
With this syntax, each "\*_R1.fastq.gz" file will be matched with the corresponding "\*_R2.fastq.gz" file.
The workflow will produce two folders and one report file: 
* the report.html file contains a report of the tasks run by the workflow with details on the resources used and outcome
* the "work" directory that contains intermediate files. These files allow the workflow manager to reuse already performed steps if something goes wrong or you want to change something (option "-resume")
* the "results" folder that contains results files. The results directory contains:
    * scaffolds: scaffolds for each input sample;
    * bins: genome bins produced by Metabat2;
    * unbinned: unbinned contigs;
    * stats_scaffolds.tsv: scaffold statistics;

See https://github.com/metashot/mag-illumina for complete documentation

### Step 2: Bin quality

```
cd quality

nextflow run metashot/prok-quality \
  --genomes '../mag/results/bins/*.fa' \
  --outdir results \
  ß-with-report report.html
```
This command uses the bin fasta files from the previous step as input, and produces a report on thei completeness and contamination. Again, it will produce a "work" directory and a "results" directory
Mail output in the "results" directory:

* genome_info.tsv: summary table of genomes quality. Columns are: 
    * Genome: the genome filename
    * Completeness, Contamination, Strain heterogeneity: CheckM estimates
    * GUNC pass: if a genome doesn't pass GUNC analysis it means it is likely to be chimeric
    * Genome size (bp), ... 
    * \# predicted genes: basic genome statistics 
    * 5S rRNA, 23S rRNA, 16S rRNA \# tRNA, # tRNA types: the number and types of rRNA and tRNA genes  respecively
* filtered genomes: a folder containg the genomes passing the qulity filter
* genome_info_filtered.tsv: same as genome_info.tsv 
* derep_info.tsv: table conating the dereplication summary. Columns are: 
    * Genome: genome filename
    * Cluster: the cluster ID (from 0 to N-1)
    * Representative: is this genome the cluster representative?
* filtered_repr: folder containing the filered and representative genomes

See https://github.com/metashot/prok-quality for complete documentation

### Step 3: Taxonomic classification of bins

This step classifies the genomes passing the qulity filter of the previous step according to the GTDB (https://gtdb.ecogenomic.org/) schema and the GTDB-Tk tookit (https://github.com/Ecogenomics/GTDBTk). 

```
cd classify 
nextflow run compmetagen/prok-classify \
  --genomes "../quality/*.fa" \
  --gtdbtk_db ./release220 \
  --outdir results
   -with-report report.html
```

Output:

* bacteria_summary.tsv: the GTDB-Tk summary for bacterial genomes;
* archaea_summary.tsv: the GTDB-Tk summary for archaeal genomes;
* bacteria_genomes: folder containing the genomes classified as bacteria by GTDB-Tk;
* archaea_genomes: folder containing the genomes classified as archaea by GTDB-Tk.

Complete documentation in https://github.com/metashot/prok-classify

#### Step 4: Functional annotation of genomes

This workflow perform gene prediction and functional annotation of prokaryotic genomes. Necessary databases are downloaded automatically

```
cd annotate
nextflow run metashot/prok-annotate \
  --genomes "../quality/results/*.fa" \
  --outdir results
```

Input: prokaryotic contig/genomes in FASTA format;
* Functional annotation using Prokka;
* Functional annotation through orthology assignment by eggNOG-mapper v2 and the eggNOG v5.0 database (optional);
* KEGG ortholog assignment by KofamScan and the KOfam database (https://www.genome.jp/tools/kofamkoala/) (optional);
* Estimates KEGG pathway completeness using Anvi'o (https://merenlab.org/software/anvio/) (optional);
