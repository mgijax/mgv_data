#
# Downloads the Alliance homology TSV file to the download directory,
# then reformats it a bit, and writes "homologies.tsv" to the output directory.
# The output file has these columns:
#	ID1	ID2	YVAL
# where ID1 and ID2 are curie IDs (eg MGI:97490 and HGNC:8620) 
# and YVAL is one the values "YN", "NY", or "YY". Example:
#
#	FB:FBgn0000017 RGD:6502032 NY
#	FB:FBgn0000017 WB:WBGene00000018 YY
#	FB:FBgn0000017 ZFIN:ZDB-GENE-020809-2 YY
#	FB:FBgn0000017 ZFIN:ZDB-GENE-100812-9 NY
#

source config.sh
source utils.sh

# download from Alliance, if updated
agrfile="${DDIR}/agrOrthology.tsv"
agrfile2="${ODIR}/homologies.json"
#
if test -e "$agrfile"
then zflag=(-z "$agrfile")
else zflag=()
fi
curl -o "$agrfile" "${zflag[@]}" "$AGR_HOM_URL"

mkdir -p "${ODIR}/homologies"
${PYTHON} getHomologies.py -d "${ODIR}/homologies" < "${agrfile}"
