import sys
import os
import re
import json
from urllib.request import urlopen

### ------------------------------------------------------------------
class UrlResolver :
    def log(self, s) :
        sys.stderr.write(s)
        sys.stderr.write("\n")

    # OVERRIDE ME.
    # rgs:
    #    gcfg - the genome's config
    #    dcfg - the data config, i.e., 
    def resolve(self, gcfg, dcfg):
        return dcfg["url"]

### ------------------------------------------------------------------
# Generates URLs for gene models and genome assembies from Ensembl.
# Examples:
# ftp://ftp.ensembl.org/pub/release-99/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz
# ftp://ftp.ensembl.org/pub/release-99/gff3/danio_rerio/Danio_rerio.GRCz11.99.chr.gff3.gz
class EnsemblResolver (UrlResolver) :
    def resolve (self, gcfg, dcfg) :
        pth = dcfg.get("remotePath", gcfg["path"])
        if dcfg["type"] == "assembly":
            dcfg["url"] = dcfg["baseUrl"] + "/release-%s/fasta/%s/dna/%s.%s.dna.toplevel.fa.gz" % \
                                      (dcfg["release"], pth, pth.capitalize(), gcfg["build"])
        elif dcfg["type"] == "models":
            dcfg["url"] = dcfg["baseUrl"] + "/release-%s/gff3/%s/%s.%s.%s.gff3.gz" % \
                                      (dcfg["release"], pth, pth.capitalize(), gcfg["build"], dcfg["release"])
        else:
            raise RuntimeError("Don't know this type: " + dcfg["type"])
        return dcfg["url"]

### ------------------------------------------------------------------
# cfg:
#    assemblyId string  GCA_000001635.9_GRCm39
# Example id -> url:
#  GCA_000001635.9_GRCm39 ->
#    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz
#  GCF_000001635.25_GRCm38.p5 ->
#    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/GCF_000001635.25_GRCm38.p5_genomic.fna.gz
class NcbiResolver (UrlResolver) :
    def resolve (self, gcfg, dcfg) :
        if dcfg["type"] == "assembly":
            ident = dcfg["assemblyId"]
            (prefix, triples, version, name) = self.parseAssemblyId(ident)
            dcfg["url"] = dcfg["baseUrl"] + ("%s/%s/%s/%s_genomic.fna.gz" % (prefix, "/".join(triples), ident, ident))
            self.log("URL: " + dcfg["url"])
            return dcfg["url"]
        else:
            raise RuntimeError("Don't know this type: " + dcfg["type"])

    # given "GCA_000001635.9_GRCm39", returns ("GCA", ["000","001","635"] ".9", "GRCm39")
    def parseAssemblyId (self, ident) :
        prefix,numeric,name = ident.split("_")
        if not prefix in ["GCA","GCF"]:
            raise RuntimeError("Invalid assembly id (bad prefix): " + ident)
        nparts = numeric.split(".")
        if len(nparts) > 2:
            raise RuntimeError("Invalid assembly id (bad numeric): " + ident)
        numeric = nparts[0]
        version = "" if len(nparts) == 1 else nparts[1]
        if len(numeric) != 9:
            raise RuntimeError("Invalid assembly id (numeric not 9 digits): " + ident)
        triples = [ numeric[i:i+3] for i in (0,3,6) ]
        return (prefix, triples, version, name)


        
### ------------------------------------------------------------------
# Gets gene models from the Alliance for non-mouse organisms
class AllianceResolver (UrlResolver) :
    SNAPSHOT_CACHE = {}
    def resolve (self, gcfg, dcfg) :
        if dcfg["type"] not in ["models","orthologs","variants"] :
            raise RuntimeError("Don't know this type:" + dcfg["type"])
        dtype = dcfg["allianceDataType"]
        dcfg["url"] = self.getUrl(gcfg, dcfg)
        return dcfg["url"]

    def getSnapshotFileList(self, gcfg, dcfg):
        rel = dcfg["release"]
        if rel in self.SNAPSHOT_CACHE:
            self.log("Reusing snapshot from cache.")
            snapshot = self.SNAPSHOT_CACHE[rel]
        else:
            self.log("Getting snapshot file at: " + dcfg["allianceSnapshotUrl"] + rel)
            fd = urlopen(dcfg["allianceSnapshotUrl"] + rel)
            snapshot = json.loads(fd.read())
            fd.close()
            self.SNAPSHOT_CACHE[rel] = snapshot
        files = snapshot["snapShot"]["dataFiles"]
        #self.log(json.dumps(files, indent=2))
        return files

    def getUrl (self, gcfg, dcfg) :
        fList = self.getSnapshotFileList(gcfg, dcfg) 
        fList = list(filter(lambda f: f["dataType"]["name"] == dcfg["allianceDataType"] and f["dataSubType"]["name"] == dcfg["provider"], fList))
        if len(fList) == 0:
            raise RuntimeError("File does not exist.")
        elif len(fList) > 1:
            raise RuntimeError("File specification is not unique.")
        f = fList[0]
        return dcfg["allianceDownloadUrl"] + f["s3Path"]

### ------------------------------------------------------------------
resolverNameMap = {
    "ncbi" : NcbiResolver,
    "ensembl" : EnsemblResolver,
    "alliance" : AllianceResolver,
}

def resolveAll (gcfg) :
    for dcfg in gcfg["tracks"]:
        cls = resolverNameMap.get( dcfg.get("source", None), UrlResolver)
        cls().resolve(gcfg, dcfg)

