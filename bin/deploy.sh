#!/usr/bin/env bash

source config.sh
source utils.sh

if [[ ${ODIR} == "" ]] ; then
    die "Please specify output in config.sh (ODIR)."
fi
if [[ ${WDIR} == "" ]] ; then
    die "Please specify web deployment directory in config.sh (WDIR)."
fi

if [[ ${ODIR} != ${WDIR} ]] ; then
    logit "Copying ${ODIR} to ${WDIR}"
    rsync -av ${ODIR} ${WDIR}
    checkExit
else
    logit "Skipping data copy as output and deployment directories are the same."
fi
logit "Copying CGI scripts."
rsync -av ./fetch.cgi ${WDIR}
checkExit
rsync -av ./fetch.py ${WDIR}
checkExit
chmod ogu+x ${WDIR}/fetch.cgi
checkExit
