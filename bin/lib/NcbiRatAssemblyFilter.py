from .AssemblyFilter import AssemblyFilter
class NcbiRatAssemblyFilter (AssemblyFilter) :
    # Looking for lines like this:
    # >NC_051348.1 Rattus norvegicus strain BN/NHsdMcwi chromosome 13, mRatBN7.2, whole genome shotgun sequence
    # >NC_051349.1 Rattus norvegicus strain BN/NHsdMcwi chromosome 14, mRatBN7.2, whole genome shotgun sequence
    # >NC_051350.1 Rattus norvegicus strain BN/NHsdMcwi chromosome 15, mRatBN7.2, whole genome shotgun sequence
    #
    # Change them to put the chromosome number up front:
    # >13 NC_051348.1 Rattus norvegicus strain BN/NHsdMcwi chromosome 13, mRatBN7.2, whole genome shotgun sequence
    def processHeaderLine(self, line) :
        if "chromosome" in line and "genomic scaffold" not in line:
            ci = line.find("chromosome") + 10
            cj = line.find(",", ci)
            c = line[ci:cj].strip()
            return ">%s %s" % (c, line[1:])
