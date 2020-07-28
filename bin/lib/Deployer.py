#
# Deployer.py
#
# Deploys the data for MGV to a web-accessible place, sets up the CGI, and .htaccess.
#

import os, stat
import sys
import json

class Deployer:
    def __init__(self, builder, type, cfg, odir, wdir, cgidir, debug=False):
        self.builder = builder
        self.log = self.builder.log
        self.type = type
        self.cfg = cfg
        self.output_rootdir = odir
        self.output_dir = os.path.join(odir, cfg["name"], type)
        self.web_rootdir = wdir
        self.web_dir = os.path.join(wdir, cfg["name"], type)
        self.cgi_dir = cgidir
        self.debug = debug
        self.builder.ensureDirectory(self.web_dir)
        self.builder.ensureDirectory(self.cgi_dir)

    def deployData (self):
        if self.output_dir == self.web_dir:
            self.log("Skipping data deployment because output and web directories are the same: " + self.output_dir)
        else:
            cmd = 'rsync -av "%s" "%s"' % (self.output_dir, os.path.dirname(self.web_dir))
            self.log("Deploying %s data with command: %s" % (self.type, cmd))
            if not self.debug:
                os.system(cmd)

    def deployCgi (self):
        # copy python script
        myDir = os.path.dirname(__file__)
        scriptname = os.path.abspath(os.path.join(myDir, '../www/fetch.py'))
        cmd = "cp -f %s %s" % (scriptname, self.cgi_dir)
        self.log("Copying python CGI script: " + cmd)
        if not self.debug:
            os.system(cmd)

        # generate CGI wrapper
        FETCH_CGI = "#!/usr/bin/env bash\n# THIS IS A GENERATED FILE. See Deployer.py\n%s %s/fetch.py --cgi --dir %s\n"
        cgi = FETCH_CGI % (sys.executable, self.cgi_dir, self.web_rootdir)
        fname = os.path.join(self.cgi_dir, "fetch.cgi")
        self.log("Generating CGI wrapper: " + fname) 
        if not self.debug:
            self.log("Opening " + fname)
            with open(fname, 'w') as ofd:
                ofd.write(cgi)
            self.log("Setting permissions on CGI.")
            os.chmod(fname, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH | stat.S_IROTH)

        # copy .htaccess file
        fname = os.path.abspath(os.path.join(myDir, '../www/apache.htaccess'))
        cmd = "cp -f %s %s/.htaccess" % (fname, self.cgi_dir)
        self.log("Copying .htaccess file: " + cmd)
        if not self.debug:
            os.system(cmd)

    def deployIndex (self):
        fnames = os.listdir(self.web_rootdir)
        subdirs = []
        for fn in fnames:
            fpath = os.path.join(self.web_rootdir, fn)
            if os.path.isdir(fpath):
                c = self.builder.getCfg(fn)
                if c and "taxonid" in c:
                    subdirs.append(fn + "/")
        subdirs.sort()
        jsubdirs = json.dumps(subdirs)
        ifn = os.path.join(self.web_rootdir, "index.json")
        self.log("Generating index file: " + ifn)
        self.log(jsubdirs)
        if not self.debug:
            with open(ifn, 'w') as ifd:
                ifd.write(jsubdirs + "\n")

    def go (self) :
        self.deployData()
        self.deployCgi()
        self.deployIndex()



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
