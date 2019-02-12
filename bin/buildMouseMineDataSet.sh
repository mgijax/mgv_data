#!/bin/bash
#
set -o pipefail
#
# ---------------------
# Echos its arguments to the log file. Prepends a datetime stamp.
#
function logit {
    if [[ ${LOGFILE} ]] ; then
      echo `date` "$*" >> ${LOGFILE}
    else
      # no log file. echo to std err
      (>&2 echo `date` "$*")
    fi
}

# ---------------------
# Logs a message and exits with error code 1.
#
function die {
    logit "$*"
    exit 1
}

# ---------------------
# If the exit code of the last command ($?) is not zero, exits with a message.
#
function checkExit {
    c=$?
    if [ $c -ne 0 ]; then
        logit "ERROR: $1" 
        exit 1
    fi  
    return 0
}

# ---------------------
function doMouse {
  logit "Mouse strains url=${mouseurl}"
  # Export GFF3 files from MouseMine for each of the sequenced strains
  python getGenomeFromMouseMine.py -d ${downloadsdir} -u "${mouseurl}"
  checkExit "Get genome from MouseMine"

  # Import genomes - chunk the files.
  for  i in $( ls ${downloadsdir} ); do
    python importGenome.py -d ${outputdir} -k ${chunksize} < ${downloadsdir}/$i
    checkExit "Import genome $i"
  done
}

# ---------------------
function doHuman {
  logit "Human url=${humanurl}"
  curl ${humanurl} | gunzip | python prepEnsembl.py -g human | python importGenome.py -x 9606 -g "H.sapiens" -k ${chunksize} -d ${outputdir}
  checkExit
}

# ---------------------
function doRat {
  logit "Rat url=${raturl}"
  curl ${raturl} | gunzip | python prepEnsembl.py -g rat | python importGenome.py -x 10116 -g "R.norvegicus" -k ${chunksize} -d ${outputdir}
  checkExit
}
# ---------------------
function usage {
  echo "Usage: bash ${scriptname} [-M MouseMineURL][-H HumanGff3FileUrl][-R RatGff3FileUrl][-d downloadsDir][-o outputDir][-k chunkSize]"
  exit -1
}

# ---------------------
scriptname="$0"
downloadsdir="./downloads"
outputdir="./output"
chunksize="4000000"
humanurl=''
raturl=''
mouseurl=''

until [ -z "$1" ]  # Until all parameters used up . . .
do
    case "$1" in
    -h)
        shift
        usage
        ;;
    -d)
        shift
        downloadsdir="$1"
        ;;
    -o)
        shift
        outputdir="$1"
        ;;
    -k)
        shift
        chunksize="$1"
        ;;
    -M)
        shift
        mouseurl="$1"
        ;;
    -H)
        shift
        humanurl="$1"
        ;;
    -R)
        shift
        raturl="$1"
        ;;
    *)
        echo "Unrecognized option:" $1
        usage
    esac
    shift
done

#
mkdir -p ${downloadsdir}
checkExit "Create downloads dir: ${downloadsdir}"

#
mkdir -p ${outputdir}
checkExit "Create output dir ${outputdir}"

#
if [[ ${humanurl} ]] ; then
  doHuman
fi

#
if [[ ${raturl} ]] ; then
  doRat
fi

#
if [[ ${mouseurl} ]] ; then
  doMouse
fi
