#?Metagenomics
#?=============================
#!What is metagenomics:
# Metagenomics is the study of genetic material recovered directly from environmental samples in order to study horizontal gene transfer as a function of environmental adaptation.

#! Questions metagenomics answers?
#Which microorganisms are present in a sample?
#What abundance (quantity) is each microorganism present in?
#How do the observed abundances change across conditions?
#What are the functions of the microorganisms that are observed?

#organisms cant survive without other organisms as such we need to treat organisms in the context of their communities.

#! Two approaches to metagenomics:
#! 1. 16S rRNA sequencing:
#Component of the prokaryotic ribosome, used due to slow rate of evolution. Its gene has both highly conserved regions as well as hypervariable regions {V1-V9 regions} which are specific to taxonomical level of a species.
#Highly economical as requires fraction of the coverage, with simple analysis methods.

#! Whose genome metagenomics
#Whole genome sequencing but for community of organisms with following caveats:
#1.Organisms may be present in wildly different abundances
#2.Sequences represent random samples 
#3.Measures are relative to total number of DNA fragments 

#! How to quantify and compare results?
#alpha diversity: is diversity within a sample site 
#beta diversity: compare diversity between sample sites 

#!Steps in an analysis:
#Reducing sequencing and PCR errors
#Processing improved sequences
#Assessing error rates
#Processing OTUs {Operational Taxonomic Unit used to classify groups of closely related individuals.}
#Calling OTUs
#Phylotypes
#Phylogenetic Analysis
#Diversity
#Alpha diversity
#Beta diversity measurements
#Population-level analysis
#Phylotype-based analysis
#Phylogeny-based analysis

#!Getting the NCBI taxonomy file 
URL=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
curl -k $URL | gunzip -c > acc2taxid.txt

#new style file contains columns : accession|accession.version|taxid|gi
#some programs still run on old style format of just : gi | taxid
#we need to convert the new style to this old style to use these programs:
cat acc2taxid.txt | grep -v taxid | awk ' { print $4,"\t", $3 } ' | head

#!How do i search the NCBI taxonomy data:
#access via wither NCBI entrez direct or taxonkit {more advances}
#whats taxid 9913 stand for?
fetch -db taxonomy -id 9913
taxonkit list --show-rank --show-name --ids 9913

#taxid of a taxonomy name 
esearch -query "E. coli[txn]" -db taxonomy | efetch
echo Escherichia coli | taxonkit name2taxid 

#can search a list of taxonomy names
cat names.txt | taxonkit name2taxid

#taxid of an accession number
esearch -db nuccore -query NC_009565,NC_002570,NC_008358 | elink -target taxonomy | efetch
efetch -db taxonomy -id 336982,272558,228405

#Check the lineage of a taxid 
efetch -db taxonomy -id 562 -format xml > out.xml
cat output.xml | xtract -pattern LineageEx/Taxon -element TaxId,ScientificName

#How many species of viruses are in current taxonomy?
taxonkit list --ids 10239 --indent "" --show-rank --show-name  | grep species | head

#! How do i get all know bacterial genomes
# Get the summary as a tabular text file.
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# Filter for complete genomes.
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths

# Identify the FASTA files (.fna.) other files may also be downloaded here.
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths

#Download the data from the urls collected:
mkdir all
cat ftpfilepaths | parallel -j 20 --verbose --progress "cd all && curl -O {}"

#to find all files in all folder with *fna.gz need use:
find all -name '*fna.gz' | wc -l

#! How to set up BLAST for taxonomy operations:
#Download taxdb 
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz

#move to the blast db directory 
blastdbcmd -db 16SMicrobial -entry all -outfmt "%a %T %t" | head

#!Extract sequences deposited under specific taxid e.g. 1392
taxonkit list --show-name --show-rank --ids 1392
URL=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
curl $URL | gunzip -c > acc2taxid.txt
# Filter accession numbers by the known taxonomy 9606.
cat acc2taxid.txt | awk ' $3 == 1392 {print $1}'  > acc.txt
# Extract sequences that belong to that taxonomy.
cat acc.txt | blastdbcmd -db nt -entry_batch -  > sequences.fa

#! Database containing enviromental samples?
mkdir -p ~/blastdb
(cd ~/blastdb && update_blastdb.pl --decompress env_nt)
blastdbcmd -db env_nt -entry all -outfmt "%a %T %S" | head

#? Classifying using 16S sequence
#? ================================
#Data: Human Microbiome Project (HMP) Demonstration Projects that contain 1065 sequencing runs that attempt to connect the disease of Necrotizing enterocolitis a devastating complication of prematurity to microbial composition of the gut.

#Get all SRR runs for the project (Total of 1065)
# Search the runs for this project.
esearch -query PRJNA46337 -db sra | efetch -format runinfo > runinfo.csv

#Download the first id SRR1614899
fastq-dump --split-files SRR1614899 -O srr

#investigate srr files metadata looking for the barcode the instrument tagged the samples with and the primer used to isolate the sequences
#In this example the sequence was TCAGACGGCTC and the 16S primer was CCGTCAATTCATTTGAGT which is the 907R primer.
#SRR1614899_3.fastq contains primary dataset used for data analysis

#! Classifying a 16S sequence
#Alignment is wasteful and redundant newer methods rely upon recognizing patterns in the sequence.
#RPD classifier is the CL choice 
java -jar ~/src/RDPTools/classifier.jar classify SRR1614899_3.fastq -o assign.txt -h ranks.txt

#the output file assign.txt has a taxonomical assignment for each input 

#! Classifying whole genome data
#Data:Assessment of Metagenomic Assembly Using Simulated Next Generation Sequencing Data, PLoS One, 2012.

#Obtain the data:
curl http://www.bork.embl.de/~mende/simulated_data/illumina_100species.1.fq.gz | gunzip -c > read1.fq
curl http://www.bork.embl.de/~mende/simulated_data/illumina_100species.2.fq.gz | gunzip -c > read2.fq

#Run analysis on read{1,2}.fq though limited it to 1000000 random reads due to size
seqtk sample -2 -s 2 read1.fq 1000000 > subset.fq

# See what the file contains.
seqkit stat subset

#looking at metadata we see Psychrobacter cryohalolentis K5 {accession numbers NC_007969} has an expected coverage of 17.84.
#Lets see if data supports this:
# Total number of reads 
cat subset.fq | egrep '^@NC_007969|^@NC_007968' | wc -l

# Get the size of NC_007969 + NC_007968
efetch -db nuccore -id NC_007969,NC_007968 -format fasta | seqkit stat
#coverage for NC_007969 = 849442 * 75 / 3101097 = 20.5x not 17.4 as stated in excel file

#estimate is 15902 * 53 ~ 842,000 total numer of reads lets check if it was close
cat read1.fq read2.fq | egrep '^@NC_007969|^@NC_007968' | wc -l
#very close at 849442

#kraken not working












