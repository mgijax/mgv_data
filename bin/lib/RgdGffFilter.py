from .AllianceGffFilter import AllianceGffFilter
class RgdGffFilter (AllianceGffFilter) :
    def processModel (self, model) :
        # Starting with Alliance 3.2.0, RGD GFF3 files (for rat and human) fixed the issues we previously had to correct for.
        # However, they also now have multiple providers of gene models (NCBI and Ensembl), but they do not merge them like we do.
        # Also, there appears to be duplication of some of the ENSEMBL models. (E.g. for PAX2).
        # As a quick solution, just pick the NCBI models. If a gene does not have an NCBI model, it gets filtered out
        # FIXME: Could/should merge the models. Or maybe get RGD to do that.

        if len(list(filter(lambda f: f[1] == "NCBI", model))) > 0:
            return AllianceGffFilter.processModel(self, model)
        else:
            return None
