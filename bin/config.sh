
# downloads directory
DDIR="./downloads"
# output directory
ODIR="./output"
# web directory (can be the same as the output directory)
WDIR="./output"
# first part of url for downloading Ensembl genomes data
ENSEMBL_BASE="rsync://ftp.ensembl.org/ensembl/pub"
# default release number to use. Can override with -r command line option.
ENSEMBL_RELEASE="97"
# Regular expression for matching chromosomes
# Default matches chromosomes of one or 2 characters
# (good for screening out scaffolds)
CHR_REGEX="..?" 
# feature types (column 3 values) to exclude
EXCLUDE_TYPES="biological_region,chromosome,scaffold"
# transcript file chunk size
CHUNK_SIZE="4000000"
# Python executable
PYTHON="python2.7"

