from .Filter import Filter
class AssemblyFilter (Filter) :
    def __init__ (self, *args):
        Filter.__init__(self, *args)
        self.showing = True

    def processNext (self, line) :
        if line.startswith(">") :
            chrom = line.split(maxsplit=1)[0][1:]
            self.showing = self.chr_re.match(chrom)
            if self.showing:
                self.log(line[:-1])
                return line
        elif self.showing:
            return line
