#!/usr/bin/env bash
#
# getGenome.sh 
#
# Downloads genome data (gene model annotations and genome assemblies) from Ensembl
# for a specified Ensembl release and organism. Imports into MGV data files. 
#
source config.sh
source utils.sh

# ---------------------
usage () {
  echo "Usage: bash ${SCRIPT} [-g organism-path][-n organism-name][-r release][-t PART][-p PHASE]"
  exit -1
}

# ---------------------
downloadModels () {
  #
  logit "Downloading ${GFF_URL} to ${GFF_GZ_FILE}."
  mkdir -p ${DDIR}
  if [[ "${GFF_URL}" == rsync* ]] ; then
      logit "Using rsync"
      rsync -av --progress ${DRY_RUN} "${GFF_URL}" "${GFF_GZ_FILE}"
  else
      logit "Using curl"
      curl -z "${GFF_GZ_FILE}" -o "${GFF_GZ_FILE}" "${GFF_URL}"
  fi
  checkExit
}

# ---------------------
importModels () {
  logit "Importing ${ORGANISM} gene models to ${ODIR}."
  mkdir -p "${G_ODIR}"
  if [[ ${DRY_RUN} == "" ]] ; then
    if [[ ("${FORCE}" != "") || ("${GFF_GZ_FILE}" -nt "${ODIR}/${ORGANISM}/index.json") ]] ; then
	gunzip -c "${GFF_GZ_FILE}" | \
	${PYTHON} prepGff3.py -x "${EXCLUDE_TYPES}" -c "${CHR_REGEX}" ${MODULES} | \
	${PYTHON} importGff3.py -p ${ORGANISM} -g ${NAME} -x ${TAXONID} -k ${CHUNK_SIZE} -d ${ODIR}
    else
        logit "Skipped import because file has not been updated. Use --force to override."
    fi
  fi
}

# ---------------------
downloadAssembly () {
  #
  logit "Downloading ${ORGANISM} genome assembly to ${FASTA_GZ_FILE}."
  mkdir -p ${DDIR}
  rsync -av --progress ${DRY_RUN} "$FASTA_URL" "${FASTA_GZ_FILE}"
  checkExit "Failed downloading ${FASTA_URL} to ${FASTA_GZ_FILE}"
}

# ---------------------
importAssembly () {
  logit "Importing ${ORGANISM} genome assembly to ${G_ODIR}/sequences."
  mkdir -p "${G_ODIR}"
  mkdir -p "${G_ODIR}/sequences"
  if [[ ${DRY_RUN} == "" ]] ; then
    if [[ ("${FORCE}" != "") || ("${GFF_GZ_FILE}" -nt "${ODIR}/${ORGANISM}/index.json") ]] ; then
	gunzip -c "${FASTA_GZ_FILE}" | ${PYTHON} importFasta.py -c "${CHR_REGEX}" -o ${G_ODIR}/sequences
    else
        logit "Skipped import because file has not been updated. Use --force to override."
    fi
  fi
}

# ---------------------
SCRIPT="$0"
MODULES=""
RELEASE="${ENSEMBL_RELEASE}"
until [ -z "$1" ]  # Until all parameters used up . . .
do
    case "$1" in
    --help|-h)
        # print help and exit
        usage
        ;;
    -g)
        # set the organism path name, eg mus_musculus_dba2j
        shift
        ORGANISM="$1"
        ;;
    -G)
        # match organism. If specified and not == ORGANISM, exit quietly.
	shift
	MATCH_ORGANISM="$1"
	;;
    -n)
        # set the organism name, eg DBA/2J
        shift
        NAME="$1"
        ;;
    --gff-url)
        # explicitly set the url to use for downloading the GFF3 file (overrides Ensembl URL)
        shift
        GFF_URL="$1"
        ;;
    --fasta-url)
        # explicitly set the url to use for downloading the Fasta file (overrides Ensembl URL)
        shift
        FASTA_URL="$1"
        ;;
    -x)
        # set the taxon id
        shift
        TAXONID="$1"
        ;;
    -r)
        # specify ensembl release number, eg 97
        shift
        RELEASE="$1"
        ;;
    -t)
        # which data type to work on (models, assembly)
        shift
        DATATYPE="$1"
        ;;
    -m)
        # set the organism path name, eg mus_musculus_dba2j
        shift
        MODULES="-m $1"
        ;;
    --exclude)
        # which SO types to exclude from GFF3 file
        shift
        EXCLUDE_TYPES="$1"
        ;;
    --chr-regex)
        # which chromsomes to process.
        shift
        CHR_REGEX="$1"
        ;;
    --force)
        # which chromsomes to process.
        FORCE="true"
        ;;
    -p)
        # which phase to run (download, import)
        shift
        PHASE="$1"
        ;;
    -D)
        # flag to specify that we're doing a dry run
        DRY_RUN="--dry-run"
        ;;
    *)
        # anything else is an error
        echo "Unrecognized option:" $1
        usage
    esac
    shift
done
#
if [[ ${DDIR} == "" ]] ; then
    die "Please specify downloads directory in config.sh (DDIR)."
fi
if [[ ${ODIR} == "" ]] ; then
    die "Please specify output directory in config.sh (ODIR)."
fi
if [[ !(${ORGANISM} && ${RELEASE}) ]] ; then
  die "Please specify both organism and release."
fi
if [[ ${MATCH_ORGANISM} && (${MATCH_ORGANISM} != ${ORGANISM}) ]] ; then
  exit 0
fi
if [[ ! ${NAME} ]] ; then
  NAME="${ORGANISM}"
fi
#
if [[ "${GFF_URL}" == "" ]] ; then
    GFF_URL="${ENSEMBL_BASE}/release-${RELEASE}/gff3/${ORGANISM}/*.${RELEASE}.gff3.gz"
fi
if [[ "${FASTA_URL}" == "" ]] ; then
    FASTA_URL="${ENSEMBL_BASE}/release-${RELEASE}/fasta/${ORGANISM}/dna/*.dna.toplevel.fa.gz"
fi
GFF_GZ_FILE="${DDIR}/${ORGANISM}.${RELEASE}.gff3.gz"
FASTA_GZ_FILE="${DDIR}/${ORGANISM}.${RELEASE}.fa.gz"
G_ODIR="${ODIR}/${ORGANISM}"
G_WDIR="${WDIR}/${ORGANISM}"
#
if [[ ${DATATYPE} == "" || ${DATATYPE} == "models" ]] ; then
  if [[ ${PHASE} == "" || ${PHASE} == "download" ]] ; then
    downloadModels
  fi
  if [[ ${PHASE} == "" || ${PHASE} == "import" ]] ; then
    importModels
  fi
fi
#
if [[ ${DATATYPE} == "" || ${DATATYPE} == "assembly" ]] ; then
  if [[ ${PHASE} == "" || ${PHASE} == "download" ]] ; then
    downloadAssembly
  fi
  if [[ ${PHASE} == "" || ${PHASE} == "import" ]] ; then
    importAssembly
  fi
fi
