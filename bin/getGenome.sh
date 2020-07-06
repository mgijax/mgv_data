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
help () {
  cat << pleh

  This script adds a genome (usually, one from Ensembl) to the data set served to MGV.
  Required command line options specify the genome name, printable label, and taxon id.
  Other options may control what part(s) of the process to run, which chromosomes to import, and other aspects.
  Adding a genome adds two data types, genome assemblies (Fasta) and genome annotations (GFF3),
  and consists of two phases, download and import.
  By default, both phases are run for both data types, but either or both can be specified
  via command line options -t and -p.

  Options:
      -g GENOME	Name of the genome used in paths, e.g., "mus_musculus_dba2j". Required.

      -G GENOME	If this option is specified and is not the same as the -g GENOME, exits quietly

      -n LABEL	Printable label for the genome, e.g., "DBA/2J"

      -x TAXON	NBCI taxon ID. Required.

      --gff-url URL	If specified, downloads the GFF3 file from here instead of Ensembl.

      --fasta-url URL	If specified, downloads the genome assembly file from here instead of Ensembl.

      -r RELEASE	If specified, overrides ENSEMBL_RELEASE specified in config.sh

      -t TYPE		One of: models, assembly. If specified, only processes this data type. By
      			default, processes both types.

      -p PHASE		One of: download, import. If specified, only runs this phase. By default,
      			runs both phases.

      --exclude-types SOTYPES	Comma-separated list of SO types to exclude doing GFF3 filtering. If specified,
      			GFF3 records with any of the specified values in the type column are filtered out.
			Overrides the values value set in config.sh.

      --chr-regex REGEX	A regular expression specifying which chromosomes to include. Applies in both 
      			GFF3 and Fasta import phases. During GFF3 import, the regex is matched against 
			the first column. During Fasta import, the regex is appied to the ID in the
			sequence header line. Overrides the value set in config.sh.

      --chunk-size SIZE	Specifies the approximate size (in bp) to use to split gene model structure data.
      			Overrides the value set in config.sh.

      -m MODULES	If specified, includes these processing filters during import of GFF3.
      			Provides a hook for doing simple transformations and filtering.
			MODULES is a comma-separated list of module names (no .py extension).
			These must exist in the bin/filter directory. Each
			defines a function named "feature" which takes a feature as argument
			and returns the feature (possibly modified), or None. If None is returned,
			the feature is omitted/skipped. The feature is a list of 9 items indexed
			from 0 to 8. The last item (the attributes column) is a dict from attr name
			to value, which is a string or list of strings. For example, here is a
			module that adds a "length" attribute to column 9:
			    # bin/filters/addLength.py
			    def feature(f) :
				f[8]["length"] = f[4] - f[3] + 1
				return f

      --force		The download and import phases use file modification times to suppress the
      			operation if the file has not changed. Specifying --force on the 
			command line will force the operation to happen regardless of modification dates.

pleh
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
        if [[ "${GFF_GZ_FILE}" =~ ^.*\.gz ]] ; then
            gunzip -c "${GFF_GZ_FILE}" | \
            ${PYTHON} prepGff3.py -x "${EXCLUDE_TYPES}" -c "${CHR_REGEX}" ${MODULES} | \
            ${PYTHON} importGff3.py -p "${ORGANISM}" -g "${NAME}" -x "${TAXONID}" -k "${CHUNK_SIZE}" -d "${ODIR}"
        else
            cat "${GFF_GZ_FILE}" | \
            ${PYTHON} prepGff3.py -x "${EXCLUDE_TYPES}" -c "${CHR_REGEX}" ${MODULES} | \
            ${PYTHON} importGff3.py -p "${ORGANISM}" -g "${NAME}" -x "${TAXONID}" -k "${CHUNK_SIZE}" -d "${ODIR}"
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
  if [[ ${DRY_RUN} == "" ]] ; then
	mkdir -p "${G_ODIR}/sequences"
	gunzip -c "${FASTA_GZ_FILE}" | ${PYTHON} importFasta.py -c "${CHR_REGEX}" -o ${G_ODIR}/sequences
	touch "${G_ODIR}/sequences"
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
        help
	exit 0 
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
    --exclude-types)
        # which SO types to exclude from GFF3 file
        shift
        EXCLUDE_TYPES="$1"
        ;;
    --chr-regex)
        # which chromsomes to process.
        shift
        CHR_REGEX="$1"
        ;;
    --chunk-size)
        # chunk size for splitting up transcript files
        shift
        CHUNK_SIZE="$1"
        ;;
    -p)
        # which phase to run (download, import)
        shift
        PHASE="$1"
        ;;
    -D)
        # flag to specify that we're doing a dry run
        DRY_RUN="true"
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
if [[ "${GFF_URL}" =~ ^.*\.gz$ ]]; then
    GFF_GZ_FILE="${DDIR}/${ORGANISM}.${RELEASE}.gff3.gz"
else
    GFF_GZ_FILE="${DDIR}/${ORGANISM}.${RELEASE}.gff3"
fi 

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
