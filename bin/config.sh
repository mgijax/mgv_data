# discover the installation directory (don't change this)
export DIR=`dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )`
# downloads directory
export DDIR="${DIR}/downloads"
# output directory
#export ODIR="${DIR}/output"
export ODIR="${DIR}/allianceOutput"
# web directory By default, is the same as the build directory.
export WDIR="${ODIR}"
# cgi script directory. By default, is the same as the web data directory.
export CDIR="${WDIR}"
# first part of url for downloading Ensembl genomes data
export ENSEMBL_BASE="rsync://ftp.ensembl.org/ensembl/pub"
# default release number to use. Can override with -r command line option.
export ENSEMBL_RELEASE="97"
# Regular expression for matching chromosomes
# Default matches chromosomes of one or 2 characters
# (good for screening out scaffolds)
# can override with --chr-regex command line option.
export CHR_REGEX="..?" 
# feature types (column 3 values) to exclude
# can override with --exclude-types command line option.
export EXCLUDE_TYPES="biological_region,chromosome,scaffold"
# transcript file chunk size
export CHUNK_SIZE="4000000"
# Python executable
export PYTHON="python2.7"
# Where to get the genome annotation file for C57BL/6J
export MGI_URL="http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz"
# Where to get the Alliance homology TSV file
export AGR_HOM_URL="http://download.alliancegenome.org/2.3.0/ORTHOLOGY-ALLIANCE/COMBINED/ORTHOLOGY-ALLIANCE_COMBINED_45.tsv"
