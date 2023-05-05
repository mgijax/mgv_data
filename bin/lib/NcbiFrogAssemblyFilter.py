from .AssemblyFilter import AssemblyFilter
class NcbiFrogAssemblyFilter (AssemblyFilter) :
    # Looking for lines like this:
    #  >CM030340.1 Xenopus laevis strain J_2021 chromosome 1L, whole genome shotgun sequence
    #  >CM030341.1 Xenopus laevis strain J_2021 chromosome 1S, whole genome shotgun sequence
    #
    # Change them to put the chromosome number up front:
    #  >1S CM030341.1 Xenopus laevis strain J_2021 chromosome 1S, whole genome shotgun sequence
    def processHeaderLine(self, line) :
        if "chromosome" in line:
            ci = line.find("chromosome") + 10
            cj = line.find(",", ci)
            c = line[ci:cj].strip()
            if 'laevis' in self.gcfg['name'] and self.gcfg['name'][-2] == ".":
                subgenome = self.gcfg['name'][-1]
                if not c.endswith(subgenome):
                    return None
            return ">%s %s" % (c, line[1:])
