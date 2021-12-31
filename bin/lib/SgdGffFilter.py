from .GffFilter import GffFilter
class SgdGffFilter (AllianceGff) : 
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

