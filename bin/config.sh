
PYTHON="python2.7"
GENOMES_FILE="genomes.tsv"

# downloads directory
DDIR="./downloads"
# output directory
ODIR="./output"
# web directory (can be the same as the output directory)
WDIR="./output"

# 
ENSEMBL_BASE="rsync://ftp.ensembl.org/ensembl/pub"

# By default, does not
MGI_MODELS="false"
MGI_URL="http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz"

ORGANISM=""
NAME=""
TAXONID=""
RELEASE=""
DRY_RUN=""
CHR_REGEX="..?"  # matches chromosomes of one or 2 characters
EXCLUDE_TYPES="biological_region,chromosome,scaffold"
K="4000000"

DATATYPE="" # models, assembly. If not specified, run all parts
PHASE="" # download, import. If not specified, run all phases

