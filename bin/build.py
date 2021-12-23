#
# build.py
#
# Builds a back end for MGV based on a config file.
#

import os
import sys
import time
import json
import yaml
from argparse import ArgumentParser
import re
from urllib.request import urlopen
import gzip

from lib.Downloader import downloaderNameMap
from lib.Importer import importerNameMap
from lib.Deployer import Deployer

### ------------------------------------------------------------------
class MgvDataBuilder :
    VALID_TYPES = ["assembly", "models", "orthologs"]
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
            "-b", "--build-config",
            required=True,
            help = "Build config file. Required.")
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
            "--snapshot-file",
            help = "Alliance release snapshot file to use in lieu of querying API. (default = get snapshot from Alliance API)")
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

    def process(self, g) :
        self.log("Processing cfg: " + str(g))
        gn = g["name"]
        for t in self.VALID_TYPES:
            if self.args.type in [t, None] :
                if not t in g:
                    continue
                #
                if type(g[t]) is str and g[t].startswith("="):
                    if "deploy" in self.args.phase:
                        gg = self.getCfg(g[t][1:])
                        tgtPath = os.path.join(self.args.web_dir, gg["name"], t)
                        lnkPath = os.path.join(self.args.web_dir, g["name"], t)
                        cmd = 'ln -s %s %s' % (tgtPath, lnkPath)
                        self.log("Creating symlink: " + cmd)
                    continue
                sname = g[t].get("source","UrlDownloader")
                cls = downloaderNameMap[sname]
                downloader = cls(self, g, t, self.args.debug)
                # Download data
                if "download" in self.args.phase:
                    downloader.go()
                # Import data
                if "import" in self.args.phase:
                    icls = importerNameMap[t]
                    importer = icls(self, t, g, self.args.output_dir, self.args.debug)
                    importer.go()
                # Deploy
                if "deploy" in self.args.phase:
                    deployer = Deployer(self, t, g, self.args.output_dir, self.args.web_dir, self.args.cgi_dir, debug=self.args.debug)
                    deployer.go()

    def getCfg (self, name = None) :
        if name is None:
            return self.cfg
        else:
            return self.name2cfg.get(name, None)

    def main (self) :
        #
        self.args = self.getArgs()
        if self.args.log_file:
            self.logfile = open(self.args.log_file, 'w')
        self.log("\n\nThis is the MGV back end data builder.")
        self.log("Arguments: " + str(self.args))
        self.genome_re = re.compile('^' + self.args.genome + '$')
        #
        self.log("Reading config file: " + self.args.build_config)
        with open(self.args.build_config, 'r') as cfile:
            self.cfg = yaml.safe_load(cfile)["buildList"]
        #self.log(json.dumps(self.cfg, indent=2))
        #
        if self.args.debug:
            self.log("Running in DEBUG mode. No commands will be executed.")
            self.log("Fully expanded build configuration:")
            self.log(yaml.safe_dump(self.cfg))
        #
        self.name2cfg = {}
        for g in self.cfg:
            self.name2cfg[g["name"]] = g
        #
        for g in self.cfg:
            if g.get("disabled", False) :
                continue
            if self.genome_re.match(g["name"]):
                self.log("Processing " + g["name"])
                self.process(g)
            else:
                # self.log("Skipping " + g["name"])
                pass
        self.log("Builder exiting.")
        self.logfile.close()


### ------------------------------------------------------------------
if __name__ == "__main__":
    MgvDataBuilder().main()
