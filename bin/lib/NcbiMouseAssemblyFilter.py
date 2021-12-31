from .AssemblyFilter import AssemblyFilter
class NcbiMouseAssemblyFilter (AssemblyFilter) :
    # Looking for lines like this:
    #    >CM000994.3 Mus musculus chromosome 1, GRCm39 reference primary assembly C57BL/6J
    # Change them to put the chromosome number up front:
    #    >1 CM000994.3 Mus musculus chromosome 1, GRCm39 reference primary assembly C57BL/6J
    # Also want MT:
    #    >AY172335.1 Mus musculus strain C57BL/6J mitochondrion, complete genome
    #
    # Don't want anything with "contig" in it:
    #    >GL456233.2 Mus musculus chromosome X unlocalized genomic contig MMCHRX_RANDOM_CTG2, GRCm39 reference primary assembly C57BL/6J
    # Or
    #    >GL456378.1 Mus musculus unplaced genomic contig MSCHRUN_CTG3, GRCm39 reference primary assembly C57BL/6J
    def processNext(self, line) :
        if line.startswith(">") and not "contig" in line:
            if "mitochondrion" in line:
                return ">MT " + line[1:]
            ci = line.find("chromosome") + 10
            cj = line.find(",", ci)
            c = line[ci:cj].strip()
            return ">%s %s" % (c, line[1:])
        return line

