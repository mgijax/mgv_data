#!/usr/bin/bash
#
# Wrapper script to launch the Python builder, which does the heavy lifting.
# 

# Discover my installation directory. Taken from:
#  https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
#
set -o pipefail

export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export INSTALL_DIR=`dirname ${SCRIPT_DIR}`

source ${SCRIPT_DIR}/wrapperConfig.sh

CMDLINE="$*"    # save for future reference
MATCHPAT=""    # genomes must match pattern to be processed (filename pattern matching is used)
GDIR=""         # path name for genome dir
PHASE="all"     # which phase to run: download, import, deploy, all
FLOCAL=""       # local name for downloaded file (= basename of the url)
TRACK="all"     # one of: models, assembly, all
FTYPE=""        # one of: gff or fasta
FURL=""         # URL of the file to import
FILTER=""       # filter program
FORCE=""        # if true, force downloads
DEBUG=""        # debug mode
LOGFILE=""      # log file

# ---------------------
usage () {
  logit "
Usage: $0 [parameters]

Parameters:
-h 
    Print this help and exit.
-g PATTERN
    Only process genomes matching the pattern. To process all genomes, -g all.
-G PATTERN
    Process all genomes, starting from the first one matching PATTERN. 
    Useful for restarting a failed build - specify the genome that failed as the argument to -G.
-t TRACK
    Which track to import: either models or assembly. Default processes all tracks.
-p PHASE
    Specifies which build phase to run. Default runs all phases. One of: download, import, deploy
-f
    If specified, force file downloads (normally, only downloads if the source has a newer modification date)..
-D
    Don't actually do anything - just print log messages of what would be done.
-L FILE
    Specifies a log file. By default, log messages are sent to stdout.
"
}

# ---------------------
parseCommandLine () {
# ---------------------
# 
    # Process command line args
    until [ -z "$1" ] 
    do
        case "$1" in
        -h) 
            # print help and exit
            usage
            exit 0
            ;;  
        -g)
            # genome path/subdir name
            shift
            MATCHPAT="$1"
            ETSEQ=""
            if [[ ${MATCHPAT} == "all" ]] ; then
                MATCHPAT="*"
            fi
            ;;
        -G)
            # genome path/subdir name
            shift
            MATCHPAT="$1"
            ETSEQ="true"
            if [[ ${MATCHPAT} == "all" ]] ; then
                MATCHPAT="*"
            fi
            ;;
        -t)
            # file type (models or assembly)
            shift
            TRACK="$1"
            ;;
        -p)
            # build phase
            shift
            PHASE="$1"
            ;;
        -f)
            # force downloads
            FORCE="true"
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
      echo `date` $D $GDIR $STAGE $FTYPE "$*" >> ${LOGFILE}
  else
      echo `date` $D $GDIR $STAGE $FTYPE "$*"
  fi
}

# ---------------------
# Logs a message and exits with error code 1.
#
die () {
    logit "$*"
    logit "Build FAILED. Exiting."
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
  if [[ ${FORCE} == "" ]] ; then
      zarg="-z ${localFile}"
  fi
  logit "${CURL} ${zarg} -o ${localFile} ${FURL}"
  if [[ $DEBUG == "" ]] ; then
    ${CURL} ${zarg} -o "${localFile}" "${FURL}"
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
  ofile="${odir}/${TTYPE}.${FTYPE}.gz"
  if [[ $FLOCAL == *.gz ]] ; then
    stream="${GUNZIP} -c"
  else
    stream="cat"
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
      if [[ $DEBUG == "" ]] ; then
          tmpFile=`mktemp ${TDIR}/mgvtmpXXXX`
          checkexit
      else
          tmpFile="TMPFILE"
      fi
      logit "$stream $ifile | ${FILTER} > ${tmpFile}"
      if [[ $DEBUG == "" ]] ; then
          $stream $ifile | ${FILTER} > ${tmpFile}
          checkexit
      fi

      logit "Compressing..."
      logit "(grep '^#' ${tmpFile} | grep -v '^###'; grep -v '^#' ${tmpFile} | sort -k1,1 -k4,4n)  | ${BGZIP} > ${ofile}"
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

      logit "Extracting top level features..."
      if [[ $DEBUG == "" ]] ; then
          gfile="${odir}/${TTYPE}.genes.gff.gz"
          ${BGZIP} -c -d ${ofile} | grep -v "^#" | grep -v "Parent=" | ${BGZIP} > ${gfile}
          ${TABIX} -p gff ${gfile}
          checkexit
      fi


  elif [ $FTYPE == "fasta" ] ; then
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

  makedirectory ${WDIR}/${GDIR}

  logit "cp -f ${ODIR}/${GDIR}/${FLOCAL}* ${WDIR}/${GDIR}"
  if [[ $ODIR == $WDIR ]] ; then
      logit "Skipped file copy because output and deployment directories are the same."
  elif [[ $DEBUG == "" ]] ; then
      cp -f ${ODIR}/${GDIR}/${TTYPE}.* ${WDIR}/${GDIR}
      checkexit
  fi

  infofile=${WDIR}/${GDIR}/index.json
  logit "Writing info file: ${infofile}"
  if [[ $DEBUG == "" ]] ; then
      echo $GCONFIG > ${infofile}
      checkexit
  fi

}

# ---------------------
deployWwwContents () {
  GDIR=""
  STAGE="deploy"
  if [[ $PHASE == "all" || $PHASE == "deploy" ]] ; then
      FTYPE="www"
      simpleCommand cp ${INSTALL_DIR}/www/apache.htaccess ${WDIR}/.htaccess
      simpleCommand cp ${INSTALL_DIR}/www/info.html ${WDIR}
      simpleCommand cp ${INSTALL_DIR}/www/fetch.py ${CDIR}
      cgi=${CDIR}/fetch.cgi
      logit "Generating: ${cgi}"
      echo "#!/usr/bin/env bash" > ${cgi}
      echo "# THIS IS A GENERATED FILE." >> ${cgi}
      echo "export TABIX=\"${TABIX}\"" >> ${cgi}
      echo "export SAMTOOLS=\"${SAMTOOLS}\"" >> ${cgi}
      echo "${PYTHON} ${CDIR}/fetch.py --cgi --dir ${WDIR} \$*" >> ${cgi}
      chmod +x ${cgi}
      logit `cat ${cgi}`
  fi
  STAGE=""
  FTYPE=""
}

# ---------------------
importGenome () {
    if [[ $FURL == "" ]] ; then
        die "No URL."
    fi

    makedirectory $TDIR

    FLOCAL=`basename $FURL`
    if [[ $PHASE == "all" || $PHASE == "download" ]] ; then
        download
    fi
    if [[ $PHASE == "all" || $PHASE == "import" ]] ; then
        import
    fi
    if [[ $PHASE == "all" || $PHASE == "deploy" ]] ; then
        deploy
    fi     
}

# ---------------------
importHomology () {
    if [[ $FURL == "" ]] ; then
        die "No URL."
    fi

    makedirectory $TDIR

    FLOCAL=`basename $FURL`
    ifile="${DDIR}/${GDIR}/${FLOCAL}"
    odir="${ODIR}/${GDIR}"
    if [[ $PHASE == "all" || $PHASE == "download" ]] ; then
        download
    fi
    if [[ $PHASE == "all" || $PHASE == "import" ]] ; then
        makedirectory ${ODIR}/${GDIR}
        logit "gunzip -c ${ifile} | ${FILTER}"
        if [[ $DEBUG == "" ]] ; then
            gunzip -c ${ifile} | ${FILTER}
            checkexit
        fi

    fi
    if [[ $PHASE == "all" || $PHASE == "deploy" ]] ; then
        FLOCAL=""
        deploy
    fi     
}


# ---------------------
activateVenv () {
    if [ -d "${SCRIPT_DIR}/venv/bin" ] ; then
        source ${SCRIPT_DIR}/venv/bin/activate
    else
        logit "Creating virtual environment."
        simpleCommand ${PYTHON} -m venv "${SCRIPT_DIR}/venv"
        source "${SCRIPT_DIR}/venv/bin/activate"
        logit "Installing yaml"
        simpleCommand pip install pyyaml
    fi
}

# ---------------------
# ---------------------
# ---------------------
