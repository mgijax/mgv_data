#!/usr/bin/env bash

source config.sh
source utils.sh

if [[ ${ODIR} == "" ]] ; then
    die "Please specify output in config.sh (ODIR)."
fi
if [[ ${WDIR} == "" ]] ; then
    die "Please specify web deployment directory in config.sh (WDIR)."
fi
logit "deploy.sh: Deploying from ${ODIR} to ${WDIR}"

if [[ ${ODIR} != ${WDIR} ]] ; then
    # rsynch the data
    logit "rsync'ing ${ODIR} to ${WDIR}"
    rsync -av ${ODIR} ${WDIR}
    checkExit
else
    logit "Skipping rsync because output and deployment directories are the same."
fi
# copy the scripts
logit "Generating CGI wrapper"
cat > ${CDIR}/fetch.cgi << EOF
#!/usr/bin/env bash
# THIS IS A GENERATED FILE. See deploy.sh
${PYTHON} ${CDIR}/fetch.py --cgi --dir ${WDIR}
EOF
checkExit
logit "Setting execute permission."
chmod ogu+x ${CDIR}/fetch.cgi
checkExit
#
rsync -av ./fetch.py ${CDIR}
checkExit
# build the root index.json file which names each of the available subdirectories.
logit "Building ${WDIR}/index.json"
cd ${WDIR}
SEP=""
echo "[" > index.json
# List subdirectories of the root that contain index.json files
for i in $(ls -F */index.json)
do
    echo $SEP "\"$(dirname ${i})/\"" >> index.json
    SEP=","
done
echo "]" >> index.json

