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

from lib.Config import ConfigFileReader
from lib.Downloader import downloaderNameMap
from lib.Importer import importerNameMap
from lib.Deployer import Deployer

### ------------------------------------------------------------------
class MgvDataBuilder :
    VALID_TYPES = ["assembly", "models", "orthology"]
    VALID_PHASES = ["download", "import", "deploy"]
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
            choices = self.VALID_PHASES,
            action = "append",
            default = [],
            help = "Which phase to run. One of: %(choices)s. If not specified, runs all phases.")
        parser.add_argument(
            "-t", "--type",
            choices = self.VALID_TYPES,
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
            "-o", "--output-dir",
            default = "./output",
            help = "Where the output files go. Default = %(default)s")
        parser.add_argument(
            "-w", "--web-dir",
            help = "Web accessible directory containing data generated files. Default = same as --output-dir.")
        parser.add_argument(
            "--cgi-dir",
            help = "Place to put the CGI scripts used by MGV Default = same as --web-dir.")
        parser.add_argument(
            "-D", "--debug",
            action = "store_true",
            default = False,
            help = "Run in debug mode.")
        args = parser.parse_args()
        args.downloads_dir = os.path.abspath(args.downloads_dir)
        args.output_dir = os.path.abspath(args.output_dir)
        args.web_dir = os.path.abspath(args.web_dir) if args.web_dir else args.output_dir
        args.cgi_dir = os.path.abspath(args.cgi_dir) if args.cgi_dir else args.web_dir
        if len(args.phase) == 0:
            args.phase = self.VALID_PHASES

        return args

    def makeDownloaderObject(self, g, type) :
        sname = g[type].get("source","UrlDownloader")
        cls = downloaderNameMap[sname]
        return cls(self, g, type, self.args.debug)

    def readConfigFile (self, fname) :
        with open(fname) as fd:
            cfg = json.load(fd)
        return cfg

    def deepCopy (self, obj) :
        return json.loads(json.dumps(obj))

    def ensureDirectory (self, d, empty = False):
        if self.args.debug:
            return
        if not os.path.exists(d):
            os.makedirs(d)
        if empty:
            cmd = "rm -fr %s/*" % d
            self.log(cmd)
            os.system(cmd)

    def process(self, gg) :
        self.log("Processing cfg: " + str(gg))
        gn = gg["name"]
        for t in self.VALID_TYPES:
            if self.args.type in [t, None] :
                if not t in gg:
                    continue
                # Download data
                downloader = self.makeDownloaderObject(gg, t)
                if "download" in self.args.phase:
                    downloader.downloadData()
                # Import data
                if "import" in self.args.phase:
                    cls = importerNameMap[t]
                    impobj = cls(self, t, gg, self.args.output_dir, self.args.debug)
                    impobj.go()
                # Deploy
                if "deploy" in self.args.phase:
                    Deployer(self, t, gg, self.args.output_dir, self.args.web_dir, self.args.cgi_dir, debug=self.args.debug).go()

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
        if self.args.debug:
            self.log("Running in DEBUG mode. No commands will be executed.")
        #
        for g in self.cfg:
            if self.genome_re.match(g["name"]):
                self.process(g)
        self.log("Builder exiting.")
        self.logfile.close()


### ------------------------------------------------------------------
if __name__ == "__main__":
    MgvDataBuilder().main()
