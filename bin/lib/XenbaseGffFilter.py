from .AllianceGffFilter import AllianceGffFilter

class XenbaseGffFilter (AllianceGffFilter) :
    def processModel (self, m) :
        # check the top-level feature in the model. Filter out models where the curie (col 9)
        # is not a Xenbase id
        f = m[0]
        attrs = f[8]
        if not attrs['curie'] or not attrs['curie'].startswith('Xenbase'):
            return None
        c = f[0]
        if 'laevis' in self.gcfg['name'] and self.gcfg['name'][-2] == ".":
            subgenome = self.gcfg['name'][-1]
            if not c.endswith(subgenome):
                return None
        return AllianceGffFilter.processModel(self, m)

    def processFeature (self, f) :
        AllianceGffFilter.processFeature(self, f)
        c = f[0].lower()
        if c.startswith('chr0'):
            f[0] = f[0][4:]
        elif c.startswith('chr'):
            f[0] = f[0][3:]
        return f
