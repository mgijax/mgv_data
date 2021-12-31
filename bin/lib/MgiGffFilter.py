from .GffFilter import GffFilter
class MgiGffFilter (GffFilter) :
    def __init__ (self, impobj) :
        GffFilter.__init__(self, impobj)

    # To allow build 38 and 39 to coexist in the viewer, need to modify the feature IDs so they're unique.
    def processModel (self, model) :
        mid = model[0][8]['ID']
        newMid = mid + '_' + GCONFIG['build']
        model[0][8]['ID'] = newMid
        for f in model[1:]:
            if mid in f[8].get('Parent', []):
                f[8]['Parent'] = [newMid]
        return GffFilter.processModel(self, model)

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

