from .Filter import Filter
class GffFilter (Filter) :
    def __init__(self, src):
        Filter.__init__(self, Gff3Parser(src, returnHeader=True, returnGroups=True).iterate())

    def processNext (self, obj) :
        if type(obj) is str:
            # the header
            return obj
        else:
            # obj is a list of features
            return ''.join(self.processModel(obj))

    def processModel(self, model):
        model = list(filter(lambda x: x, [self._processFeature(f) for f in model]))
        '''
        #universal transform: make sure top level feature's coordinates span its descendants
        if len(model) > 1:
            f = model[0] # top level feature
            f[3] = min([c[3] for c in model[1:]])
            f[4] = max([c[4] for c in model[1:]])
        '''
        return model

    def _processFeature (self, f):
        # universal transform: filter out features on non-matching chromosome 
        if not CHR_RE.match(f[0]):
            return None
        # universal transform: filter out features with excluded types
        if 'include_types' in DCONFIG and f[2] not in DCONFIG['include_types']:
            return None
        if 'exclude_types' in DCONFIG and f[2] in DCONFIG['exclude_types']:
            return None
        return formatLine(self.processFeature(f))

