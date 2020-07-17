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
import gzip

from lib.Downloader import downloaderNameMap
from lib.Importer import GffImporter, FastaImporter
from lib.Config import ConfigFileReader

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

    def makeDownloaderObject(self, g, type) :
        sname = g[type].get("source","UrlDownloader")
        cls = downloaderNameMap[sname]
        return cls(self, g, type)

    def readConfigFile (self, fname) :
        with open(fname) as fd:
            cfg = json.load(fd)
        return cfg

    def deepCopy (self, obj) :
        return json.loads(json.dumps(obj))

    def processGenome (self, gg) :
        gn = gg["name"]
        for t in ["assembly", "models"] :
            if self.args.type in [t, None] :
                downloader = self.makeDownloaderObject(gg, t)
                # Download genome data
                if self.args.phase in ["download", None] :
                    self.log("%s: downloading %s: %s" % (gn, t, downloader.cfg[t]["url"]))
                    downloader.downloadData()
                # Import genome data
                if self.args.phase in ["import", None] :
                    if t == "assembly":
                        impobj = FastaImporter(self, "assembly", gg, self.args.output_dir)
                    elif t == "models":
                        impobj = GffImporter(self, "models", gg, self.args.output_dir)
                    self.log("%s: importing %s" % (gn, t))
                    impobj.go()

    def main (self) :
        #
        self.args = self.getArgs()
        if self.args.log_file:
            self.logfile = open(self.args.log_file, 'w')
        self.log("\n\nThis is the MGV back end data builder.")
        self.log("Arguments: " + str(self.args))
        self.genome_re = re.compile('^' + self.args.genome + '$')
        #
        self.cfg = ConfigFileReader(self.args.config_file).read()
        for g in self.cfg:
            if self.genome_re.match(g["name"]):
                self.log(str(g))
                self.processGenome(g)
            else:
                self.log("Skipped %s." % g["name"])
        self.log("Builder exiting.")
        self.logfile.close()


### ------------------------------------------------------------------
if __name__ == "__main__":
    MgvDataBuilder().main()
