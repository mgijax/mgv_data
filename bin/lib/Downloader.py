import os
import re
import json
from urllib.request import urlopen

### ------------------------------------------------------------------
class Downloader :
    def __init__ (self, builder, cfg, type) :
        self.cfg = cfg
        tcfg = cfg[type]
        self.builder = builder
        self.log = self.builder.log
        self.type = type
        self.init()
        if "url" not in tcfg:
            raise RuntimeError("%s: url is not set." % self.__class__.__name__)
        # Download directory. Each genome gets its own subdirectory.
        self.cfg["ddir"] = os.path.abspath(os.path.join(self.builder.args.downloads_dir, self.cfg["name"]))
        # Set the local downloaded file name to be the same as in the URL
        fname = tcfg["url"].split("/")[-1]
        tcfg["fpath"] = os.path.abspath(os.path.join(self.cfg["ddir"], fname))
        #
        if "chr_re" in self.cfg:
            self.cfg["chr_re"] = re.compile(self.cfg["chr_re"])
        else:
            self.cfg["chr_re"] = re.compile(".*")

    def runCommand (self, cmd) :
        self.log("Running command: " + cmd)
        os.system(cmd)

    def init(self):
        pass

    def downloadData (self) :
        tcfg = self.cfg[self.type]
        if not "url" in tcfg:
            raise RuntimeError("No URL")
        #
        if tcfg["url"].startswith("rsync://"):
           cmd = 'rsync -av --progress "%s" "%s"' % (tcfg["url"], tcfg["fpath"])
        else:
           cmd = 'curl -z "%s" -o "%s" "%s"' % (tcfg["fpath"], tcfg["fpath"], tcfg["url"])
        if not self.builder.args.debug:
            self.runCommand("mkdir -p %s" % self.cfg["ddir"])
            self.runCommand(cmd)

### ------------------------------------------------------------------
class UrlDownloader (Downloader) :
    pass

### ------------------------------------------------------------------
# Gets gene models and genome assembies from Ensembl.
# Examples:
# ftp://ftp.ensembl.org/pub/release-99/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz
# ftp://ftp.ensembl.org/pub/release-99/gff3/danio_rerio/Danio_rerio.GRCz11.99.chr.gff3.gz
class EnsemblDownloader (Downloader) :
    BASEURL = "rsync://ftp.ensembl.org/ensembl/pub"
    def init (self) :
        c = self.cfg
        t = self.cfg[self.type]
        pth = c.get("remotePath", c["name"])
        if self.type == "assembly":
            t["url"] = self.BASEURL + "/release-%s/fasta/%s/dna/%s.%s.dna.toplevel.fa.gz" % \
                                      (t["release"], pth, pth.capitalize(), c["build"])
        elif self.type == "models":
            t["url"] = self.BASEURL + "/release-%s/gff3/%s/%s.%s.%s.gff3.gz" % \
                                      (t["release"], pth, pth.capitalize(), c["build"], t["release"])
        else:
            raise RuntimeError("Don't know this type: " + self.type)

### ------------------------------------------------------------------
# Gets gene models from the Alliance for non-mouse organisms
class AllianceDownloader (Downloader) :
    SNAPSHOT_CACHE = {}
    SNAPSHOT_URL="https://fms.alliancegenome.org/api/snapshot/release/" 
    DOWNLOAD_URL="https://download.alliancegenome.org/"
    def init (self) :
        if self.type not in ["models","orthology"] :
            raise RuntimeError("Don't know this type:" + self.type)
        dtype = self.cfg[self.type]["dataType"]
        self.cfg[self.type]["url"] = self.getUrl(dtype)
        self.log("URL: " + self.cfg[self.type]["url"])

    def getSnapshotFileList(self):
        rel = self.cfg[self.type]["release"]
        if rel in self.SNAPSHOT_CACHE:
            self.log("Reusing snapshot from cache.")
            snapshot = self.SNAPSHOT_CACHE[rel]
        else:
            self.log("Getting snapshot file at: " + self.SNAPSHOT_URL+rel)
            fd = urlopen(self.SNAPSHOT_URL+rel)
            snapshot = json.loads(fd.read())
            fd.close()
            self.SNAPSHOT_CACHE[rel] = snapshot
        files = snapshot["snapShot"]["dataFiles"]
        return files

    def getUrl (self, type) :
        provider = self.cfg[self.type]["provider"]
        fList = self.getSnapshotFileList()
        fList = list(filter(lambda f: f["dataType"]["name"] == type and f["dataSubType"]["name"] == provider, fList))
        if len(fList) == 0:
            raise RuntimeError("File does not exist.")
        elif len(fList) > 1:
            raise RuntimeError("File specification is not unique.")
        f = fList[0]
        return self.DOWNLOAD_URL + f["s3Path"]

### ------------------------------------------------------------------
# Gets gene models for C57BL/6J mouse strain
class MgiDownloader (Downloader) :
    def init (self) :
        if self.type != "models" :
            raise RuntimeError("Don't know this type:" + self.type)
        self.cfg[self.type]["url"] = "http://www.informatics.jax.org/downloads/mgigff3/MGI.gff3.gz"

### ------------------------------------------------------------------
downloaderNameMap = {
    "ensembl" : EnsemblDownloader,
    "alliance" : AllianceDownloader,
    "mgi" : MgiDownloader,
    "urlDownloader" : UrlDownloader
}

