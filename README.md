## Elixir Course, Bari, Italy, 2025

### Hands-on n.1 - Taxonomic and functional profiling using shotgun data
#### Topic n.1: Preprocessing
##### Step n.1: get into the right place
```
cd /home/user<YOUR USER NAME>
```
##### Step n.2: set up anaconda and check whether your environment is visible
```
DON'T INSTALL IT...
##wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
##bash Anaconda3-2024.10-1-Linux-x86_64.sh

path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
source ${path}/activate

conda info --envs
```
##### Step n.3: raw data pre-processing on fastq example files "seq_1.fastq.gz" and "seq_2.fastq.gz" from https://github.com/biobakery/biobakery/wiki/kneaddata
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

##### Step n.4: Define variable "s" with the sampleID and run TRIMMOMATIC
```
s="seq"

source ${path}/activate trimmomatic

trimmomatic PE -threads 8 -phred33 -trimlog ${s}_trimmomatic.log ${s}1.fastq ${s}2.fastq \
${s}_filtered_1.fastq ${s}_unpaired_1.fastq ${s}_filtered_2.fastq ${s}_unpaired_2.fastq \
ILLUMINACLIP:${path}/../envs/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75

for i in *.fastq; do echo -ne "${i}\t"; cat "$i" | wc -l; done
```

##### Step n. 5: Generate bowtie2 index of the human genome GCF_009914755.1_T2T-CHM13v2.0.fna (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0.fna)
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
##### Did the preprocessing produce the same exact number of reads in R1 and R2 ?

#### Topic n.2: MetaPhlAn 4: taxonomic profiling using marker genes
##### Step n.1: Setup correct variables, activate environment and navigate to the right folders
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

##### Step n.2: Download metagenomic samples
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

##### Step n.3: Let's have a look at the MetaPhlAn parameters
```
metaphlan -h
```

##### Step n.4: Run MetaPhlAn 4
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

##### Step n.6:
##### Step n.7:


#### Topic n.3: Kraken + Bracken: taxonomic profiling using k-mers
#### Topic n.4: HUMAnN 4: functional profiling at the community level




### Hands-on n.2 - Taxonomic profiling beyond the level of species
#### Step n.1:

#### Step n.2: Getting example files (6 fastq files) from https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS013951.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS014613.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS019161.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS022137.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS055982.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS064276.fastq.bz2
```

#### Step n.3: Running MetaPhlAn 4
##### Approach n. 1 ==> Running MetaPhlAn 4 to obtain the .sam files of the marker genes' alignments
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

##### Approach n. 2 ==> Copy the .sam alignments from a pre-existing repository
```
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS013951.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS014613.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS019161.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS022137.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS055982.sam.bz2 .
cp /home/ubuntu/course_backup/course/4_strainphlan/SRS064276.sam.bz2 .
```

#### Step n.4: Extract for each sample the alignments over its markers
```
mpa_database="/home/ubuntu/shotgun_course/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl"
sample2markers.py -i *.sam.bz2 -o ./ -n 8 -d ${mpa_database}
```

#### Step n.5: Extract marker genes for a species of interest
```
mkdir -p db_markers
```
##### Approach n. 1 ==> run the dedicate command
```
extract_markers.py -c t__SGB1877 -o db_markers/ -d ${mpa_database} ## TOO LONG,
```

##### Approach n. 2 ==> Copy the pre-built marker files
```
cp /home/ubuntu/course_backup/course/4_strainphlan/db_markers/t__SGB1877.fna db_markers/
```

#### Step n.6: Also include a reference genome ("GCF000273725")
```
mkdir -p reference_genomes
wget -P reference_genomes/ http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/reference_genomes/G000273725.fna.bz2
```

#### Step n.7: Let's look the StrainPhlAn params
```
strainphlan -h
```

#### Step n.8: Run StrainPhlAn 4
```
mkdir -p strainphlan_output
strainphlan -s *.json.bz2 -m db_markers/t__SGB1877.fna -r reference_genomes/G000273725.fna.bz2 -o strainphlan_output -c t__SGB1877 -n 8 -d ${mpa_database}
```

#### Step n.9: Let's visualize it ! 
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/metadata.txt
add_metadata_tree.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre -f metadata.txt -m subjectID --string_to_remove .fastq.bz2

conda deactivate
source ${path}/activate graphlan
${path}/../envs/mpa/bin/plot_tree_graphlan.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre.metadata -m subjectID
```



### Hands-on n.3 - Metagenome assembly and binning
#### Step n.1:
#### Step n.2:
#### Step n.3:
#### Step n.4:
