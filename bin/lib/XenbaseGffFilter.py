import sys
from .GffFilter import GffFilter

class XenbaseGffFilter (GffFilter) :
    def processFeature (self, f) :
        attrs = f[8]
        # chromosome
        c = f[0].lower()
        if c.startswith('chr0'):
            f[0] = f[0][4:]
        elif c.startswith('chr'):
            f[0] = f[0][3:]
        # feature type
        gbt = attrs.get('gene_biotype','')
        if gbt == 'protein_coding':
            f[2] = 'protein_coding_gene'
        elif gbt == 'pseudogene':
            f[2] = 'pseudogene'
        elif gbt.endswith('RNA'):
            f[2] = 'ncRNA_gene'
        elif gbt.endswith('_segment'):
            f[2] = 'gene_segment'

        # curie, for top level features
        if not 'Parent' in attrs:
            xbid = list(filter(lambda s:s.startswith('Xenbase:'), attrs.get('Dbxref',[])))
            if len(xbid) == 1:
                attrs['curie'] = xbid[0]

        # deal with laevis S and L subgenomes
        c = f[0]
        if 'laevis' in self.gcfg['name'] and self.gcfg['name'][-2] == ".":
            subgenome = self.gcfg['name'][-1]
            if not c.endswith(subgenome):
                return None

        # May 2023: The xenopus laevis model file contains 23 CDS features for gene Xenbase:XB-GENE-17337693 where the start..end 
        # coordinates are 1733597287..86540353.
        # Filter any CDS where end < start.
        if f[2] == "CDS" and f[4] < f[3]:
            sys.stderr.write("Skipped: " + str(f) + "\n")
            return None

        #
        return f
