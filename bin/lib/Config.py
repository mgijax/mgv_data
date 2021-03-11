import sys
import json
from datetime import date

class ConfigFileReader (dict) :
    def __init__(self, fname):
        self.cfgFileName = fname

    def read (self) :
        with open(self.cfgFileName) as fd:
            cfg = json.load(fd)
        return self.expand(cfg)

    def expand(self, cfg) :
        self.vars = cfg.get("vars",{})
        # inject special var with today's date
        self.vars['@today'] = date.today().strftime('%d %b %Y')
        #
        self._expand(self.vars)
        self.defaults = cfg.get("defaults", {})
        self._expand(self.defaults)
        self.data = cfg.get("data", [])
        self.expanded = []
        for g in self.data:
            gg = self.deepCopy(self.defaults)
            gg.update(g)
            self._expand(gg)
            self.expanded.append(gg)
        return self.expanded

    def _expand (self, o) :
        if type(o) is type({}):
            # expand each dictionary item
            for k,v in list(o.items()):
                if type(v) is str and v.startswith("@"):
                    o[k] = self.deepCopy(self.vars[v[1:]])
                else:
                    self._expand(v)
            # if an item has key == "@", copy its items into the parent obj
            if '@' in o and type(o['@']) is type({}):
                base = o.pop('@', {})
                for k,v in list(base.items()):
                    if k not in o:
                        o[k] = v
        elif type(o) is type([]):
            for i,v in enumerate(o):
                if type(v) is str and v.startswith("@"):
                    o[i] = self.deepCopy(self.vars[v[1:]])
                else:
                    self._expand(v)

    def deepCopy (self, obj) :
        return json.loads(json.dumps(obj))

if __name__ == "__main__":
    print(json.dumps(ConfigFileReader(sys.argv[1]).read(), indent=2))
