#
# Filter.py
#
# Herein lies all the ugliness needed to convert GFF3 files as downloaded from their
# respective sources into GFF3 as needed for import. 
# A filter implements a set of transformations on a gene model.
# In genomes.json, each stanza specifies (by name) zero or more filters to apply during import.
# The names used in genome.json are defined in filterNameMap (at the end of this file), which maps
# each name to the class that implements it.
# A filter class is a subclass of GffFilter that defines one or both methods
# processModel and processFeature. 
# Each model is a list of individual features (gene, transcripts, exons, CDSs, etc)
# that make up the model. The general flow is that processModel is called with the
# list of features, and it then calls processFeature on each list element in turn.
# The top level feature (the gene or pseudogene) is always the first
# element of the list. 
#
# Overriding processFeature:
# - good for simple line-by-line transforms
# 

import sys
import os
import re
from .gff3lite import parseLine, formatLine
from urllib.request import urlopen

MOUSEMINE_URL=os.environ["MOUSEMINE_URL"]
CURIE_INFO = [{
  "prefix" : "ENSEMBL",
  "baseRegex" : re.compile(r'ENS[A-Z]+[GTP]+\d+'),
}, {
  "prefix" : "ENSEMBL",
  "baseRegex" : re.compile(r'MGP_[A-Z0-9]+_[GTP]\d+'),
}, {
  "prefix" : "RefSeq",
  "baseRegex" : re.compile(r'[NX][RMP]_\d+')
}]
def curie_ize (ident) :
    i = ident.find(":")
    if i >= 0:
        return ident
    for info in CURIE_INFO:
        if info["baseRegex"].match(ident):
            return "%s:%s" % (info["prefix"], ident)
    return None


class Filter:
    def __init__ (self, impobj):
        self.importer = impobj
        self.log = self.importer.log

    # OVERRIDE me
    def processObj(self, obj):
        return obj

    def __call__(self, src):
        for obj in src:
            obj = self.processObj(obj)
            if obj:
                yield obj

class GffFilter (Filter) :
    def processObj (self, obj) :
        if type(obj) is str:
            # the header
            return obj
        else:
            # obj is a list of features
            return self.processModel(obj)

    def processModel(self, model):
        model = list(filter(lambda x: x, [self._processFeature(f) for f in model]))
        #universal transform: make sure top level feature's coordinates span its descendants
        if len(model) > 1:
            f = model[0] # top level feature
            f[3] = min([c[3] for c in model[1:]])
            f[4] = max([c[4] for c in model[1:]])
        return model

    def _processFeature (self, f):
        # universal transform: filter out features on non-matching chromosome 
        if not self.importer.chr_re.match(f[0]):
            return None
        # universal transform: filter out features with excluded types
        if f[2] in self.importer.exclude_types:
            return None
        return self.processFeature(f)

class EnsemblMouseFilter (GffFilter) :
    # index mapping ensembl IDs to MGI ids
    # Shared by all instances. The first one to try to access the index creates it.
    EID2MGI = None
    def getEid2MgiIndex (self) :
        if self.EID2MGI:
            return self.EID2MGI
        self.EID2MGI = {}
        # URL for a query that returns all mouse genes with their Ensembl IDs. Data returned 
        # as TSV with 3 columns: MGIid, symbol, Ensembl id
        url = MOUSEMINE_URL + '''/mousemine/service/query/results?query=%3Cquery+name%3D%22%22+model%3D%22genomic%22+view%3D%22Gene.primaryIdentifier+Gene.symbol+Gene.crossReferences.identifier%22+longDescription%3D%22%22+sortOrder%3D%22Gene.crossReferences.identifier+asc%22+constraintLogic%3D%22A+and+B%22%3E%3Cconstraint+path%3D%22Gene.crossReferences.source.name%22+code%3D%22A%22+op%3D%22%3D%22+value%3D%22Ensembl+Gene+Model%22%2F%3E%3Cconstraint+path%3D%22Gene.dataSets.name%22+code%3D%22B%22+op%3D%22%3D%22+value%3D%22Mouse+Gene+Catalog+from+MGI%22%2F%3E%3C%2Fquery%3E&format=tab''' 

        self.log("Getting MGI/Ensembl ID associations from: " + url)
        for r in urlopen(url):
            rr = r.decode('utf-8').strip().split()
            self.EID2MGI[rr[2]] = rr
        return self.EID2MGI

    def stripPrefix(self, eid) :
        parts = eid.split(':', 1)
        if parts[0] in ["gene","transcript","CDS"]:
            return parts[1]
        return eid

    def processFeature(self, feat):
        attrs = feat[8]
        if feat[2] == "gene" and attrs.get('biotype', None) == 'protein_coding':
            feat[2] = 'protein_coding_gene'
        if 'ID' in attrs:
            attrs['ID'] = self.stripPrefix(attrs['ID'])
        if 'Parent' in attrs:
            attrs['Parent'] = list(map(self.stripPrefix, attrs['Parent']))
        if "projection_parent_gene" in attrs:
            eid = attrs['projection_parent_gene'].split(".")[0]
            mgi = self.getEid2MgiIndex().get(eid, None)
            if mgi:
                attrs['cID'] = mgi[0]
                attrs['Name'] = mgi[1]
            else:
                attrs['cID'] = eid
                attrs['Name'] = eid
        if "transcript_id" in attrs:
            attrs["transcript_id"] = "ENSEMBL:" + attrs["transcript_id"]
        if "protein_id" in attrs:
            attrs["protein_id"] = "ENSEMBL:" + attrs["protein_id"]
        return feat

class EnsemblNonMouseFilter (GffFilter) :
    # used for finding/extracting canonical id from another attribute, by taxon id
    taxon2re = {
        "7955":  ('description', re.compile(r'Source:ZFIN.*Acc:([A-Z0-9-]+)'), "ZFIN:"),
        "9606":  ('description', re.compile(r'Acc:HGNC:(\d+)'),                "HGNC:"),
        "10116": ('description', re.compile(r'Source:RGD.*Acc:(\d+)'),         "RGD:"),
        "559292":('description', re.compile(r'Source:SGD.*Acc:(S\d+)'),        "SGD:"),
        "7227":  ('gene_id',     re.compile(r'(.*)'),                          "FB:"), 
        "6239":  ('gene_id',     re.compile(r'(.*)'),                          "WB:"),

    }
    def __init__(self, *args):
        GffFilter.__init__(self, *args)
        self.cfg = self.importer.cfg
        self.attrName, self.id_re, self.id_prefix = self.taxon2re[self.cfg["taxonid"]]

    def processFeature (self, f) :
        attrs = f[8]
        if f[2] == "gene" and attrs.get('biotype', None) == 'protein_coding':
            f[2] = 'protein_coding_gene'
        attrVal = attrs.get(self.attrName,'')
        m = self.id_re.search(attrVal)
        if m:
            attrs['cID'] = self.id_prefix + m.group(1)
        return f
        
class MgiGff (GffFilter) :
    def __init__ (self, impobj) :
        Filter.__init__(self, impobj)
        self.cfg = self.importer.cfg

    # To allow build 38 and 39 to coexist in the viewer, need to modify the feature IDs so they're unique.
    def processModel (self, model) :
        mid = model[0][8]['ID']
        newMid = mid + '_' + self.cfg['build']
        model[0][8]['ID'] = newMid
        for f in model[1:]:
            if mid in f[8].get('Parent', []):
                f[8]['Parent'] = [newMid]
        return GffFilter.processModel(self, model)

    def processFeature(self, f):
        attrs = f[8]
        if f[2] in ["gene","pseudogene"]:
            f[2] = attrs['so_term_name']
            attrs['cID'] = attrs['curie']
            attrs['long_name'] = attrs.pop('description')
        if 'transcript_id' in attrs:
            attrs['transcript_id'] = curie_ize(attrs['transcript_id'])
        if 'protein_id' in attrs:
            attrs['protein_id'] = curie_ize(attrs['protein_id'])
        return f

class AllianceGff (GffFilter) : 
    def processFeature(self, f) :
        attrs = f[8]
        attrs.pop("description", None)
        if "curie" in attrs and "Parent" not in attrs:
            attrs["cID"] = attrs["curie"]
        if "so_term_name" in attrs:
            f[2] = attrs.pop("so_term_name")
        elif "Ontology_term" in attrs:
            soid = None
            soids = list(filter(lambda i: i.startswith("SO:"), attrs["Ontology_term"]))
            if len(soids):
                soid = soids[0]
            if soid == "SO:0001217":
                f[2] = "protein_coding_gene"
            elif soid == "SO:0000336":
                f[2] = "pseudogene"
            elif soid == "SO:0001263":
                f[2] = "ncRNA_gene"
        return f

#
class RgdGff (AllianceGff) :
    def processModel (self, model) :
        # Starting with Alliance 3.2.0, RGD GFF3 files (for rat and human) fixed the issues we previously had to correct for.
        # However, they also now have multiple providers of gene models (NCBI and Ensembl), but they do not merge them like we do.
        # Also, there appears to be duplication of some of the ENSEMBL models. (E.g. for PAX2).
        # As a quick solution, just pick the NCBI models. If a gene does not have an NCBI model, it gets filtered out
        # FIXME: Could/should merge the models. Or maybe get RGD to do that.

        if len(list(filter(lambda f: f[1] == "NCBI", model))) > 0:
            return AllianceGff.processModel(self, model)
        else:
            return None
#
class ZfinGff (AllianceGff) :
    def processFeature (self, f) :
        AllianceGff.processFeature(self, f)
        attrs = f[8]
        if "curie" in attrs and "Parent" in attrs:
            attrs["transcript_id"] = attrs.pop("curie")
        if "full_name" in attrs:
            attrs["long_name"] = attrs.pop("full_name")
        if f[2] == "CDS" and attrs['ID'].startswith('CDS:'):
            attrs['protein_id'] = curie_ize(attrs['ID'][4:])
        return f

#
class FlyBaseGff (AllianceGff) :
    def processFeature (self, f) :
        # Actually nothing to do. FB transcripts already have both a curie and transcript_id.
        # The GFF does not provide protein ids.
        return AllianceGff.processFeature(self, f)

#
class WormBaseGff (AllianceGff) :
    def processFeature (self, f) :
        AllianceGff.processFeature(self, f)
        attrs = f[8]
        if 'ID' in attrs and attrs['ID'].startswith('Transcript:'):
            attrs['transcript_id'] = 'WB:' + attrs['ID'][11:]
        return f

#
class SgdGff (AllianceGff) : 
    # SGF gff issues:
    # Yeast transcription/translation is simpler b.c. no introns. Just an mRNA and a CDS.
    # Therefore in the GFF:
    # - CDS features do not have IDs (no need to tie multiple CDSs together)
    #   => assign IDs based on Parent's ID
    # - CDS features can preceed their mRNA parents (forward reference issue)
    #   => have to re-sort the model's features
    # Other things:
    # - a gene's symbol (if it exists) is in the "gene" attribute
    # - chromosomes beging with "chr" in the GFF but not in the Fasta.
    #
    def processModel(self, model) :
        # re-sort by: level, then start pos. 
        def keyfun(f):
            l1 = 0
            if "Parent" in f[8]:
                l1 = 1 if f[2] in ["transcript", "mRNA"] else 2
            return (l1, f[4])
        model.sort(key=keyfun)
        exons = []
        self.hasCDS = False
        for f in model:
            if f[2] in ["transcript","mRNA"]:
               e = f[:]
               e[2] = "exon"
               e[8] = { "Parent" : f[8]["ID"] }
               exons.append(e)
            if f[2] == "CDS":
                self.hasCDS = True
                pid = f[8]["Parent"][0]
                cid = pid.replace("mRNA", "CDS")
                if cid == pid:
                    cid = pid + "_CDS"
                f[8]["ID"] = cid
        model = model + exons
        return AllianceGff.processModel(self, model)

    def processFeature(self, f):
        AllianceGff.processFeature(self, f)
        attrs = f[8]
        if f[0].startswith("chr") :
            f[0] = f[0][3:]
            if f[0] == "mt":
                f[0] = "Mito"
        if f[2] == "gene":
            if self.hasCDS:
                f[2] = "protein_coding_gene"
            if "gene" in attrs:
                attrs["Name"] = attrs.pop("gene")
                attrs["long_name"] = attrs.pop("display")
        return f

class NcbiMouseAssemblyFilter (Filter) :
    # Looking for lines like this:
    #    >CM000994.3 Mus musculus chromosome 1, GRCm39 reference primary assembly C57BL/6J
    # Change them to put the chromosome number up front:
    #    >1 CM000994.3 Mus musculus chromosome 1, GRCm39 reference primary assembly C57BL/6J
    # Also want MT:
    #    >AY172335.1 Mus musculus strain C57BL/6J mitochondrion, complete genome
    #
    # Don't want anything with "contig" in it:
    #    >GL456233.2 Mus musculus chromosome X unlocalized genomic contig MMCHRX_RANDOM_CTG2, GRCm39 reference primary assembly C57BL/6J
    # Or
    #    >GL456378.1 Mus musculus unplaced genomic contig MSCHRUN_CTG3, GRCm39 reference primary assembly C57BL/6J
    def processObj(self, line) :
        if line.startswith(">") and not "contig" in line:
            if "mitochondrion" in line:
                return ">MT " + line[1:]
            ci = line.find("chromosome") + 10
            cj = line.find(",", ci)
            c = line[ci:cj].strip()
            return ">%s %s" % (c, line[1:])
        return line

# This maps filter names used in the config.json file to classes that implement the filters
filterNameMap = {
  "ensemblMouseFilter" : EnsemblMouseFilter,
  "ensemblNonMouseFilter" : EnsemblNonMouseFilter,
  "mgiGff" : MgiGff,
  "allianceGff" : AllianceGff,
  "sgdGff" : SgdGff,
  "rgdGff" : RgdGff,
  "zfinGff" : ZfinGff,
  "flybaseGff" : FlyBaseGff,
  "wormbaseGff" : WormBaseGff,
  "ncbiMouseAssemblyFilter" : NcbiMouseAssemblyFilter,
}
