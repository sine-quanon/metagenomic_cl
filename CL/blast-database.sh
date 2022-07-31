# Error handling.
set -uex

# The directory storing the blast databases.
mkdir -p ~/db

# This should work but does not for us. Fails at a decompression step.
# (cd ~/db && update_blastdb.pl Betacoronavirus --decompress)

# Download the latest blast database with wget
(cd ~/db && parallel -j 1 curl --silent -O https://ftp.ncbi.nlm.nih.gov/blast/db/Betacoronavirus.{}.tar.gz ::: 00 01 02 03)

# Unpack the files.
(cd ~/db && ls -1 *.gz | parallel tar zxvf {})

# Download the latest representative genomes into the blast directory.
(cd ~/db && update_blastdb.pl ref_viruses_rep_genomes --decompress)

# Download the taxonomical database.
(cd ~/db && update_blastdb.pl taxdb --decompress)

# Make the blast database universally available.
export BLASTDB=~/db

# View the content of the blast database.
blastdbcmd -info -db Betacoronavirus

# How many complete genomes in the Blast database?
blastdbcmd -db Betacoronavirus -entry all -outfmt "%t" | grep "complete genome" | wc -l
