#!/bin/bash

source utils.sh

ODIR=$1
DDIR=$2

if [[ ${ODIR} == "" || $DDIR == "" ]] ; then
    die "usage: $0 output-dir deploy-dir"
fi
if [[ ${ODIR} != ${DDIR} ]] ; then
    rsync -av ${ODIR} ${DDIR}
fi
rsync -av ./fetch.cgi ${DDIR}
chmod ogu+x ${DDIR}/fetch.cgi
rsync -av ./sequenceHound.py ${DDIR}
