from .AllianceGffFilter import AllianceGffFilter
from .GffFilter import curie_ize
class ZfinGffFilter (AllianceGffFilter) :
    def processFeature (self, f) :
        AllianceGffFilter.processFeature(self, f)
        attrs = f[8]
        if "curie" in attrs and "Parent" in attrs:
            attrs["transcript_id"] = attrs.pop("curie")
        if "full_name" in attrs:
            attrs["long_name"] = attrs.pop("full_name")
        if f[2] == "CDS" and attrs['ID'].startswith('CDS:'):
            attrs['protein_id'] = curie_ize(attrs['ID'][4:])
        return f
