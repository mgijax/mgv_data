#
# Deployer.py
#
# Deploys the data for MGV to a web-accessible place, sets up the CGI, and .htaccess.
#

import os, stat
import sys
import json
import time

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

    # Returns the date of the most recent update to any file in the given directory tree
    # Date is returned as a string.
    def getMostRecentUpdate (self, d) :
        maxt = 0
        for directory, subdirs, files in os.walk(d):
            for f in files:
                ff = os.path.join(directory, f)
                s = os.stat(ff)
                mt = s.st_mtime
                maxt = max(maxt, mt)
        mtt = time.asctime(time.localtime(maxt))
        return mtt

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

    def getChromsomesAndLengthsFromModels (self) :
        # scan the genes file. For each chroosome, keep the max end coordinate
        fpath = os.path.join(self.web_rootdir, self.cfg['name'], 'models', 'genes', '0.gff3')
        if not os.path.isfile(fpath):
            return []
        chroms = {}
        with open(fpath, 'r') as fd:
            for line in fd:
                fields = line.split('\t')
                chrom = fields[0]
                coord = int(fields[4])
                chroms[chrom] = max(chroms.get(chrom,0), coord)
        chromsList = []
        for chrom, length in chroms.items() :
            chromsList.append({ "name":chrom, "length": length })
        return chromsList

    def getChromsomesAndLengthsFromAssembly (self):
        directory = os.path.join(self.web_rootdir, self.cfg['name'], 'assembly')
        if not os.path.isdir(directory):
            return []
        chrs = []
        for f in os.listdir(directory):
            fpath = os.path.join(directory, f)
            if os.path.isfile(fpath) :
                c = f.split('.')[0]
                s = os.stat(fpath)
                sz = s.st_size
                chrs.append({ "name" : c, "length" : sz })
        return chrs

    def getChromsomesAndLengths (self) :
        c_a = self.getChromsomesAndLengthsFromAssembly()
        c_b = self.getChromsomesAndLengthsFromModels()
        n2ca = dict(map(lambda ca: (ca["name"], ca), c_a))
        n2cb = dict(map(lambda cb: (cb["name"], cb), c_b))
        for cb in c_b:
            if not cb["name"] in n2ca:
                c_a.append(cb)
                n2ca[cb["name"]] = cb
            else:
                ca = n2ca[cb["name"]]
                ca["length"] = max(ca["length"], cb["length"])
        chrs = c_a
        csort = self.cfg.get("chr_sort","standard")
        chrs.sort(key = romanSortKey if csort == "roman" else standardSortKey)
        return chrs

    def deployGenomeIndexFile (self) : 
        timestamp = self.getMostRecentUpdate(self.output_dir)
        self.log(timestamp, self.output_dir)
        ixfile = os.path.join(self.web_rootdir, self.cfg['name'], 'index.json')
        self.log("Writing genome info to " + ixfile)
        c = self.cfg
        self.log(str(c))
        tcs = c["models"]["transcriptChunkSize"]
        info = {
          "name" : c["label"],
          "timestamp" : timestamp,
          "chromosomes" : self.getChromsomesAndLengths(),
          "tracks": [
            {   
              "name": "genes",
              "type": "ChunkedGff3",
              "chunkSize": 0
            },  
            {   
              "name": "transcripts",
              "type": "ChunkedGff3",
              "chunkSize": tcs
            },  
            {   
              "name": "sequences",
              "type": "PlainSequence"
            }   
          ],  
          "linkouts" : c["models"].get("linkouts",[]),
          "metadata" : {
              "taxonid" : c["taxonid"],
              "assemblyBuild" : c["build"],
              "annotationSource" : c["models"]["source"],
              "annotationRelease" : c["models"]["release"],
          }
        }
        with open(ixfile, 'w') as fd:
            s = json.dumps(info, indent = 2)
            fd.write(s)


    def go (self) :
        self.deployData()
        if 'taxonid' in self. cfg:
            self.deployGenomeIndexFile()
        self.deployIndex()
        self.deployCgi()


def standardSortKey (c) :
    c = c['name'].upper()
    if c.isdigit() :
        return (int(c), '')
    else:
        return (9999999, c)

def romanSortKey (c) :
    c = c['name'].upper()
    n = parseRoman(c)
    if n:
        return (n, '')
    else:
        return (9999999, c)

def parseRoman (s):
      roman = {'I':1,'V':5,'X':10,'L':50,'C':100,'D':500,'M':1000,'IV':4,'IX':9,'XL':40,'XC':90,'CD':400,'CM':900}
      i = 0
      num = 0
      try:
          while i < len(s):
             if i+1<len(s) and s[i:i+2] in roman:
                num+=roman[s[i:i+2]]
                i+=2
             else:
                #print(i)
                num+=roman[s[i]]
                i+=1
          return num
      except:
          return None

