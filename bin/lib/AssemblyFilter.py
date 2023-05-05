from .Filter import Filter
class AssemblyFilter (Filter) :
    def __init__ (self, *args):
        Filter.__init__(self, *args)
        self.showing = True

    # A Fasta sequence header line (beginning with ">") has been detected.
    # Subclasses can override to apply custom modifications to the line.
    # Return None to skip the sequence.
    def processHeaderLine (self, line):
        return line

    # Process next Fasta line. If it's header line, check that the identifier after the ">" matches
    # the configured regex. If not, skip the sequence. If so, call the header processor method (which
    # may modify the line). If it's happy, output the header line and the rest of the sequence.
    def processNext (self, line) :
        if line.startswith(">") :
            chrom = line.split(maxsplit=1)[0][1:]
            self.showing = self.chr_re.match(chrom)
            if self.showing:
                line = self.processHeaderLine(line)
                if line:
                    self.log(line[:-1])
                    return line
                else:
                    self.showing = False
        elif self.showing:
            return line
