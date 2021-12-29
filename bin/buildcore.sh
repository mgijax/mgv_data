#!/usr/bin/bash
#
# Wrapper script to launch the Python builder, which does the heavy lifting.
# 

# Discover my installation directory. Taken from:
#  https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
#
export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source ${SCRIPT_DIR}/wrapperConfig.sh

# ---------------------
# 
MATCHPAT="*"
CMDLINE="$*"
PHASE=""
DEBUG=""
LOGFILE=""
FLOCAL=""
FTYPE=""
FURL=""
FILTER="cat"
GDIR=""

TIMESTAMP=`date`
STAGE=""

# ---------------------
function usage {
  logit "
Usage: $0 

Parameters:
-h 
    Print this help and exit.
-g GENOME
    Required. Genome path/output directory (e.g. mus_musculus_grcm39), under ODIR (specified in wrapperConfig.sh).
-m PATTERN
    Only do anything if GENOME matches PATTERN.
-t TYPE
    Required. Specifies the data type being imported, either gff or fa.
-p PHASE
    Optional. Specifies which build phase to run. Default runs all phases. One of: download, import, deploy
-u URL
    Required. Specifies the external URL for the GFF or FA file. May be compressed or uncompressed. 
    For local files, use URL beginning with file://
-l NAME
    Optional. Specify a local name for the file. Default is the basename from the URL.
-f COMMAND
    Optional. Filter command. Before import, file is piped through COMMAND to apply any
    needed transformations/filters.
-D
    Don't actually do anything - just print log messages of what would be done.
-L FILE
    Specifies a log file. By default, log messages are sent to stdout.
"
}

# ---------------------
function parseCommandLine {
    # Process command line args
    until [ -z "$1" ] 
    do
        case "$1" in
        -h) 
            # print help and exit
            usage
            exit 0
            ;;  
        -u)
            # specify resource URL
            shift
            FURL="$1"
            ;;
        -l)
            # specify local file name (optional; default is the basename from $FURL
            shift
            FLOCAL="$1"
            ;;
        -g)
            # genome path/subdir name
            shift
            GDIR="$1"
            ;;
        -m)
            # genome path/subdir name
            shift
            MATCHPAT="$1"
            ;;
        -t)
            # file type (gff or fa)
            shift
            FTYPE="$1"
            ;;
        -p)
            # build phase
            shift
            PHASE="$1"
            ;;
        -f)
            # filter command 
            shift
            FILTER="$1"
            ;;
        -D)
            # debug mode. Just logging (no changes).
            DEBUG="true"
            ;;
        -L)
            # logfile. (default writes to stdout)
            shift
            LOGFILE="$1"
            ;;
        *)
            # Error
            die "Unrecognized option: $1"
        esac
        shift
    done

    if [[ $FURL == "" ]] ; then
        die "No URL. Please specify -u."
    fi
    if [[ $FLOCAL == "" ]] ; then
        FLOCAL=`basename $FURL`
    fi
}

# ---------------------
# Echos its arguments to the log file. Prepends a datetime stamp.
#
logit () {
  D="(D)"
  if [[ $DEBUG == "" ]] ; then
    D=""
  fi
  if [[ ${LOGFILE} ]] ; then
      echo `date` $D $STAGE "$*" >> ${LOGFILE}
  else
      echo `date` $D $STAGE "$*"
  fi
}

# ---------------------
# Logs a message and exits with error code 1.
#
die () {
    logit "$*"
    exit 1
}

# ---------------------
# If the exit code of the last command ($?) is not zero, exits with a message.
#
checkexit () {
    c=$?
    if [ $c -ne 0 ]; then
        if [ "$1" == "-w" ] ; then
            shift
            logit "WARNING: nonzero exit code detected." "$*"
        else
            die "ERROR: $*"
        fi
    fi
    return 0
}

# ---------------------
# Logs a command, performs it, and checks the error code.
# E.g.:
#        simpleCommand mkdir -p foo/bar
# Or
#        python prog.py -a foo -b bar
simpleCommand () {
    logit "$*"
    if [[ $DEBUG == "" ]] ; then
        $*
        checkexit
    fi
}

# ---------------------
makedirectory () {
    logit "mkdir -p $1"
    if [[ $DEBUG == "" ]] ; then
        mkdir -p $1
        checkexit
    fi
}
# ---------------------
download () {
  STAGE="download"
  #
  #makedirectory ${DDIR}/${GDIR}
  simpleCommand mkdir -p ${DDIR}/${GDIR}

  localFile="${DDIR}/${GDIR}/${FLOCAL}"
  logit "${CURL} -z ${localFile} -o ${localFile} ${FURL}"
  if [[ $DEBUG == "" ]] ; then
    ${CURL} -z "${localFile}" -o "${localFile}" "${FURL}"
    checkexit
  fi
}

# ---------------------
import () {
  STAGE="import"
  #
  ifile="${DDIR}/${GDIR}/${FLOCAL}"
  #
  odir="${ODIR}/${GDIR}"
  if [[ $FLOCAL == *.gz ]] ; then
    stream="${GUNZIP} -c"
    ofile="${odir}/${FLOCAL}"
  else
    stream="cat"
    ofile="${odir}/${FLOCAL}".gz
  fi
  #
  if [[ ${FILTER} == "" ]] ; then
      FILTER="cat"
  fi
  #
  makedirectory ${ODIR}/${GDIR}
  #
  if [ $FTYPE == "gff" ] ; then
      logit "Filtering..."
      tmpFile=`mktemp ${TDIR}/mgvtmpXXXX`
      checkexit
      logit "$stream $ifile | ${FILTER} > ${tmpFile}"
      if [[ $DEBUG == "" ]] ; then
          $stream $ifile | ${FILTER} > ${tmpFile}
          checkexit
      fi

      logit "Compressing..."
      logit "(grep '^#' ${tmpFile} ; grep -v '^#' ${tmpFile} | sort -k1,1 -k4,4n)  | ${BGZIP} > ${ofile}"
      if [[ $DEBUG == "" ]] ; then
          (grep "^#" ${tmpFile} ; grep -v "^#" ${tmpFile} | sort -k1,1 -k4,4n)  | ${BGZIP} > ${ofile}
          checkexit
      fi
      #
      logit "rm -f ${tmpFile}"
      rm -f ${tmpFile}
      checkexit

      logit "Creating tabix index..."
      logit "${TABIX} -p gff ${ofile}"
      if [[ $DEBUG == "" ]] ; then
          ${TABIX} -p gff ${ofile}
          checkexit
      fi
  elif [ $FTYPE == "fa" ] ; then
      logit "Filtering and compressing..."
      logit "$stream $ifile | ${FILTER} | ${BGZIP} > ${ofile}"
      if [[ $DEBUG == "" ]] ; then
          $stream $ifile | ${FILTER} | ${BGZIP} > ${ofile}
          checkexit
      fi

      logit "Creating faidx index..."
      logit "${FAIDX} ${ofile}"
      if [[ $DEBUG == "" ]] ; then
          ${FAIDX} ${ofile}
          checkexit
      fi
  fi
}

# ---------------------
deploy () {
  STAGE="deploy"
  logit "Deploying..."
  if [[ $ODIR == $WDIR ]] ; then
      logit "Skipped file copy because output and deployment directories are the same."
      return
  fi

  makedirectory ${WDIR}/${GDIR}

  logit "cp -f ${ODIR}/${GDIR}/${FLOCAL}* ${WDIR}/${GDIR}"
  if [[ $DEBUG == "" ]] ; then
      cp -f ${ODIR}/${GDIR}/${FLOCAL}* ${WDIR}/${GDIR}
      checkexit
  fi
}

# ---------------------
function main {
    STAGE="start"
    #
    parseCommandLine $*
    if [[  $GDIR != $MATCHPAT ]] ; then
        return
    fi
    if [[ $PHASE == "" || $PHASE == "download" ]] ; then
        download
    fi
    if [[ $PHASE == "" || $PHASE == "import" ]] ; then
        import
    fi
    if [[ $PHASE == "" || $PHASE == "deploy" ]] ; then
        deploy
    fi     
}

# ---------------------
# ---------------------
# ---------------------
