#!/bin/bash
#
# Build script for preparing data files for MGV.
#
# Example of a full build:
# bash build.sh \
#   -M http://www.mousemine.org/mousemine/service \
#   -o ./output \
#   -d ./downloads
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
  # for each available genome in MouseMine
  for  gname in $( python getGenomeFromMouseMine.py ); do
    gfname=`echo "${gname}" | tr '[:upper:]' '[:lower:]' | tr -d ' /'`
    gfpath="${downloadsdir}/${gfname}.gff3"
    # Export GFF3 file for this genome
    python getGenomeFromMouseMine.py -g "${gname}" -d - -u "${mouseurl}" > "${gfpath}"
    checkExit "Get genome ${gname} from MouseMine"
    #
    python importGff3.py -d ${outputdir} -k ${chunksize} -c "##sequence-region" < ${gfpath}
    checkExit "Import genome ${gname}"
  done
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
mouseurl=""

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
if [[ ${mouseurl} ]] ; then
  doMouse
fi
