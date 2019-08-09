source config.sh

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
