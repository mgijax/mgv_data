#
# Deployer.py
#
# Deploys the data for MGV to a web-accessible place, sets up the CGI, and .htaccess.
#

import os
import sys

class Deployer:
    def __init__(self, builder, type, cfg, odir, wdir):
        self.builder = builder
        self.log = self.builder.log
        self.type = type
        self.cfg = cfg
        self.output_dir = odir
        self.web_dir = wdir
        self.cgi_dir = web_dir

    def deployData (self):
        if self.output_dir == self.web_dir:
            self.log("Skipping data deployment because output and web directories are the same: " + self.output_dir)
        else:
            cmd = 'rsync -av "%s" "%s"' % (self.output_dir, self.web_dir)
            self.log("Deploying data with command: " + cmd)

    def deployCgi (self):
        self.log("Deploying CGI wrapper")
        cgi = FETCH_CGI % (sys.executable, self.web_dir, self.data_dir)

    def deployIndex (self):
        pass

    def go (self) :
        self.deployData()
        self.deployCgi()
        self.deployIndex()

FETCH_CGI = '''#!/usr/bin/env bash
# THIS IS A GENERATED FILE. See Deployer.py
%s %s/fetch.py --cgi --dir %s
'''


'''
#!/usr/bin/env bash
#
# deploy.sh
#
# Deploys data, CGIs, and .htaccess to web-accessible locations.
#
source config.sh
source utils.sh

###################################
if [[ ${ODIR} == "" ]] ; then
    die "Please specify output in config.sh (ODIR)."
fi
if [[ ${WDIR} == "" ]] ; then
    die "Please specify web deployment directory in config.sh (WDIR)."
fi
###################################

logit "deploy.sh: Deploying from ${ODIR} to ${WDIR}"

###################################
# Copy the output data files to the web directory (if they are not the same)
if [[ ${ODIR} != ${WDIR} ]] ; then
    # rsynch the data
    logit "rsync'ing ${ODIR} to ${WDIR}"
    rsync -av ${ODIR} ${WDIR}
    checkExit
else
    logit "Skipping rsync because output and deployment directories are the same."
fi

###################################
# CGI scripts
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
###################################
# index.json
# Build the root file which names each of the available subdirectories.
logit "Building ${WDIR}/index.json"
pushd ${WDIR}
SEP=""
echo "[" > index.json
# List subdirectories of the root that contain index.json files
for i in $(ls -F */index.json)
do
    echo $SEP "\"$(dirname ${i})/\"" >> index.json
    SEP=","
done
echo "]" >> index.json
popd

###################################
# .htaccess
rsync -av ./apache.htaccess "${ODIR}/.htaccess"

###################################
# success
logit "deploy.sh: finished."
exit 0
'''
