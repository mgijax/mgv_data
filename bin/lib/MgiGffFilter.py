from .GffFilter import GffFilter, curie_ize
class MgiGffFilter (GffFilter) :
    def processFeature(self, f):
        attrs = f[8]
        if f[2] in ["gene","pseudogene"]:
            f[2] = attrs['so_term_name']
            attrs['long_name'] = attrs.pop('description')
        if 'transcript_id' in attrs:
            attrs['transcript_id'] = curie_ize(attrs['transcript_id'])
        if 'protein_id' in attrs:
            attrs['protein_id'] = curie_ize(attrs['protein_id'])
        return f

