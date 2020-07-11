#
# build.py
#
# Builds a back end for MGV based on a config file.
#

import os
import sys
import time
import json
from argparse import ArgumentParser
import re
from urllib.request import urlopen

### ------------------------------------------------------------------
class Source :
    def __init__ (self, args) :
        for (k,v) in args.items() :
            setattr(self, k, v)
        self.init()
        #
        self.ddir = os.path.abspath(os.path.join(self.builder.args.downloads_dir, self.name))
        fname = self.url.split("/")[-1]
        self.fpath = os.path.abspath(os.path.join(self.ddir, fname))

    def runCommand (self, cmd) :
        self.builder.log("Running command: " + cmd)
        os.system(cmd)

    def init(self):
        pass

    def downloadData (self) :
        if not self.url:
            raise RuntimeError("No URL")
        #
        if self.url.startswith("rsync://"):
           cmd = 'rsync -av --progress "%s" "%s"' % (self.url, self.fpath)
        else:
           cmd = 'curl -z "%s" -o "%s" "%s"' % (self.fpath, self.fpath, self.url)
        if not self.debug:
            self.runCommand("mkdir -p %s" % self.ddir)
            self.runCommand(cmd)

    def importData (self) :
        self.builder.log("Importing file: " + self.fpath)

### ------------------------------------------------------------------
class UrlSource (Source) :
    pass

### ------------------------------------------------------------------
# Gets gene models and genome assembies from Ensembl.
# Examples:
# ftp://ftp.ensembl.org/pub/release-99/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz
# ftp://ftp.ensembl.org/pub/release-99/gff3/danio_rerio/Danio_rerio.GRCz11.99.chr.gff3.gz
class EnsemblSource (Source) :
    BASEURL = "rsync://ftp.ensembl.org/ensembl/pub"
    def init (self) :
        if self.type == "assembly":
            self.url = self.BASEURL + "/release-%s/fasta/%s/dna/%s.%s.dna.toplevel.fa.gz" % \
                                      (self.release, self.name, self.name.capitalize(), self.build)
        elif self.type == "models":
            self.url = self.BASEURL + "/release-%s/gff3/%s/%s.%s.%s.chr.gff3.gz" % \
                                      (self.release, self.name, self.name.capitalize(), self.build, self.release)
        else:
            raise RuntimeError("Don't know this type: " + self.type)

### ------------------------------------------------------------------
# Gets gene models from the Alliance for non-mouse organisms
class AllianceSource (Source) :
    SNAPSHOT_CACHE = {}
    SNAPSHOT_URL="https://fms.alliancegenome.org/api/snapshot/release/" 
    DOWNLOAD_URL="https://download.alliancegenome.org/"
    def init (self) :
        if self.type != "models" :
            raise RuntimeError("Don't know this type:" + self.type)
        self.url = self.getUrl("GFF")

    def getSnapshotFileList(self):
        if self.release in self.SNAPSHOT_CACHE:
            snapshot = self.SNAPSHOT_CACHE[self.release]
        else:
            fd = urlopen(self.SNAPSHOT_URL+self.release)
            snapshot = json.loads(fd.read())
            fd.close()
            self.SNAPSHOT_CACHE[self.release] = snapshot
        files = snapshot["snapShot"]["dataFiles"]
        return files

    def getUrl (self, type) :
        fList = self.getSnapshotFileList()
        fList = list(filter(lambda f: f["dataType"]["name"] == type and f["dataSubType"]["name"] == self.provider, fList))
        if len(fList) == 0:
            raise RuntimeError("File does not exist.")
        elif len(fList) > 1:
            raise RuntimeError("File specification is not unique.")
        f = fList[0]
        return self.DOWNLOAD_URL + f["s3Path"]

### ------------------------------------------------------------------
# Gets gene models for C57BL/6J mouse strain
class MgiSource (Source) :
    def init (self) :
        if self.type != "models" :
            raise RuntimeError("Don't know this type:" + self.type)
        self.url = "http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz"

### ------------------------------------------------------------------
sourceNameMap = {
    "ensembl" : EnsemblSource,
    "alliance" : AllianceSource,
    "mgi" : MgiSource
}

### ------------------------------------------------------------------
class MgvDataBuilder :
    def __init__ (self) :
        self.logfile = sys.stderr
        self.genome_re = None

    def log(self, s, newline='\n', timestamp=True) :
        if timestamp:
            ts = time.asctime(time.localtime(time.time()))
            self.logfile.write(ts + " ")
        self.logfile.write(str(s))
        self.logfile.write(newline)
        self.logfile.flush()

    def getArgs (self) :
        parser = ArgumentParser("Builds the backend for MGV based on a config file.")
        parser.add_argument(
            "-c", "--config-file",
            default = "./genomes.json",
            help = "Build config file. Default = %(default)s.")
        parser.add_argument(
            "-g", "--genome",
            default = ".*",
            help = "Which genomes to build. By default, builds all genomes. Specify a regex pattern used to match the genome names.")
        parser.add_argument(
            "-p", "--phase",
            choices = ["download", "import"],
            default = None,
            help = "Which phase to run. One of: %(choices)s. If not specified, runs all phases.")
        parser.add_argument(
            "-t", "--type",
            choices = ["models", "assembly"],
            default = None,
            help = "Which datatype to process. One of: %(choices)s. If not specified, processes all types.")
        parser.add_argument(
            "-l", "--log-file",
            default = None,
            help = "Where to write log messages. By default, logs to stderr.")
        parser.add_argument(
            "-d", "--downloads-dir",
            default = "./downloads",
            help = "Where downloaded files go. Default = %(default)s")
        parser.add_argument(
            "-w", "--working-dir",
            default = "./work",
            help = "Where working files go during a build. Default = %(default)s")
        parser.add_argument(
            "-o", "--output-dir",
            default = "./output",
            help = "Where the output files go. Default = %(default)s")
        parser.add_argument(
            "-D", "--debug",
            action = "store_true",
            default = False,
            help = "Run in debug mode.")
        args = parser.parse_args()
        args.downloads_dir = os.path.abspath(args.downloads_dir)
        args.working_dir = os.path.abspath(args.working_dir)
        args.output_dir = os.path.abspath(args.output_dir)

        return args

    def readConfigFile (self, fname) :
        with open(fname) as fd:
            cfg = json.load(fd)
        return cfg

    def makeSourceObject(self, type, g) :
        if "source" in g[type]:
            sname = g[type]["source"]
            cls = sourceNameMap[sname]
        else:
            sname = "UrlSource"
            cls = Source
        # creates args from defaults (if any) plus everything in g
        args = self.cfg["defaults"].get(sname, {})
        args.update(g)
        args.update(g[type])
        args["type"] = type
        args["builder"] = self
        args["debug"] = self.args.debug
        return cls(args)

    def processGenome (self, g) :
        gn = g["name"]
        for t in ["assembly", "models"] :
            if self.args.type in [t, None] :
                srcobj = self.makeSourceObject(t, g)
                if self.args.phase in ["download", None] :
                    self.log("%s: downloading %s: %s" % (gn, t, srcobj.url))
                    srcobj.downloadData()
                if self.args.phase in ["import", None] :
                    self.log("%s: importing %s" % (gn, t))
                    srcobj.importData()

    def main (self) :
        #
        self.args = self.getArgs()
        if self.args.log_file:
            self.logfile = open(self.args.log_file, 'w')
        self.log("\n\nThis is the MGV back end data builder.")
        self.log("Arguments: " + str(self.args))
        self.genome_re = re.compile(self.args.genome)
        self.log(str(self.args))
        #
        self.cfg = self.readConfigFile(self.args.config_file)
        for g in self.cfg['genomes']:
            if self.genome_re.match(g["name"]):
                self.processGenome(g)
            else:
                self.log("Skipped %s." % g["name"])
        self.logfile.close()


### ------------------------------------------------------------------
if __name__ == "__main__":
    MgvDataBuilder().main()
