#
# getDataFromMouseMine.py
#
# Generates files for the Multiple Genome Viewer from the strain data in MouseMine.
#
# Usage:
#    % python getDataFromMouseMine.py OUTPUT_DIR
#
import types
import sys
import urllib
import os.path
import json
import time
import itertools
import argparse

# Web service address root for mousemine
MOUSEMINE="http://www.mousemine.org/mousemine/service"
# Maximum length of a feature. Features longer than this are filtered out and reported.
MAX_SIZE=10000000
#
TIMESTAMP=time.asctime(time.localtime(time.time()))

def getOpts () :
  parser = argparse.ArgumentParser(description='Generate files for MGV from MouseMine.')
  parser.add_argument('-u', '--url', dest='url', 
                    default=MOUSEMINE,
                    help='Base web service URL. (default=%(default)s)')
  parser.add_argument('-o', '--output-dir', dest='odir', 
                    default='./output',
                    help='The output directory. (default=%(default)s)')
  parser.add_argument('-m', '--max-feature-size', dest='maxSize', 
                    type=int,
                    default=MAX_SIZE,
                    help='The maximum allowed feature size. (default=%(default)d)')
  parser.add_argument('--sample', dest='sample', 
                    default=False,
                    action='store_true',
                    help='Generate sample output (default=%(default)s)')
  return parser.parse_args()

class DataGetter :
    def __init__ (self, opts) :
        self.odir = opts.odir
        self.url = opts.url
        self.max_size = opts.maxSize
        self.sample = opts.sample
        #
	#
        self.lFd = sys.stderr

    def log(self, s):
        self.lFd.write(s)
        self.lFd.write('\n')

    def doQuery (self, q) :
        fmt = 'tab'
        url = '%s/query/results?format=%s&query=%s' % (self.url,fmt,urllib.quote_plus(q))
        fd = urllib.urlopen(url)
        for line in fd:
            toks = line[:-1].split('\t')
            yield toks
        fd.close()
        
    def convertStrainName (self, s):
        if "PAHARI" in s:
            return "mus_pahari"
        elif "CAROLI" in s:
            return "mus_caroli"
        elif "SPRET" in s:
            return "mus_spretus"
        else:
            return "mus_musculus_" + s.lower().replace('/','')
          
    # Returns iterator over the set of annotated genomes.
    # Args:
    #    none
    # Returns:
    #    List of rows, one per genome.
    #    Each row has a primary identifier, a name, and a filename prefix.
    def getGenomes (self) :
        q = '''<query 
        model="genomic"
        view="Strain.primaryIdentifier Strain.name"
        sortOrder="Strain.name ASC"
        >
          <constraint path="Strain" op="IN" value="Annotated strains"/>
        </query>'''

        for r in self.doQuery(q):
            r.append(self.convertStrainName(r[1]))
            yield { 'name': r[1], 'filename': self.convertStrainName(r[1]) }

    # Returns iterator over chromosomes for the given genome
    # Args:
    #    g (string) the name of the genome, eg, "A/J"
    # Returns:
    #    List of rows, one per chromosome of the given genome.
    #    Each row has a chromosome id and total chromosome length.
    #    The rows are sorted by chromosome id.
    def getChromosomes (self, g) :
        q = '''<query 
        model="genomic"
        view="Chromosome.primaryIdentifier Chromosome.length"
        sortOrder="Chromosome.symbol ASC"
        >
        <constraint path="Chromosome.strain.name" op="=" value="%s"/>
        </query>''' % g

        for r in self.doQuery(q):
            yield {
              'name': r[0],
              'length': int(r[1])
            }

    # Returns an iterator over genes on the specified chromosome for the specified genome.
    # Args:
    #    g (string) the name of the genome, eg, "A/J"
    #    c (string) the chromosome, eg, "13"
    # Returns:
    #    Iterator that yields one row per gene in the specified genome.
    #    Each row has: strain name, gene id, chromosome location (chr, start, end, strand)
    def getGenes (self, g, c) :
        q = '''<query 
        model="genomic"
        view="
          Gene.strain.name
          Gene.primaryIdentifier
          Gene.sequenceOntologyTerm.name
          Gene.chromosome.primaryIdentifier
          Gene.chromosomeLocation.start
          Gene.chromosomeLocation.end
          Gene.chromosomeLocation.strand
          Gene.canonical.primaryIdentifier
          Gene.canonical.symbol
          "
        sortOrder="Gene.chromosomeLocation.start ASC"
        >
        <join path="Gene.canonical" style="OUTER"/>
        <constraint path="Gene.strain.name" op="=" value="%s"/>
        <constraint path="Gene.chromosome.primaryIdentifier" op="=" value="%s"/>
        </query>''' % (g, c)

        for r in self.doQuery(q):
            rr = {
              'ID': r[1],
              'sotype': r[2],
              'chr': r[3],
              'start': int(r[4]),
              'end': int(r[5]),
              'strand': '-' if r[6] == '-1' else '+',
              'cID': '.' if r[7] == '""' else r[7],
              'symbol': '.' if r[8] == '""' else r[8]
            }
            if rr['end'] - rr['start'] > self.max_size:
                self.log('Feature too big (skipped): ' + str(rr))
                continue
            yield rr

    # Returns an iterator over models for genes on the specified chromosome for the specified genome.
    # Args:
    #    g (string) the name of the genome, eg, "A/J"
    #    c (string) the chromosome, eg, "13"
    # Returns:
    #   Iterator yielding one object per gene, with each object containing that gene's model (transcripts and exons)
    def getModels(self, g, c) :
        q = '''
          <query
            model="genomic"
            view="Exon.transcripts.gene.primaryIdentifier Exon.transcripts.primaryIdentifier Exon.chromosomeLocation.start Exon.chromosomeLocation.end"
            sortOrder="Exon.transcripts.gene.primaryIdentifier asc Exon.transcripts.primaryIdentifier asc Exon.chromosomeLocation.start asc"
            constraintLogic="A and B">
              <constraint path="Exon.strain.name" code="A" op="=" value="%s"/>
              <constraint path="Exon.transcripts.gene.chromosome.primaryIdentifier" code="B" op="=" value="%s"/>
          </query>
        ''' % (g, c)
        for gid, recs in itertools.groupby(self.doQuery(q), lambda r: r[0]):
          rr = {
            'gID' : gid,
            'transcripts' : []
          }
          for tid, exons in itertools.groupby(recs, lambda r: r[1]):
            tr = {
              'tID' : tid,
              'exons' : [ { 'start':int(e[2]), 'end':int(e[3]) } for e in exons ]
            }
            rr['transcripts'].append(tr)
          yield rr

    def ensureDir(self, dir):
       if not os.access(dir, os.R_OK|os.W_OK):
           os.makedirs(dir)

    # Main program. 
    def main (self):
	#
	self.ensureDir(self.odir)
	# Timestamp file
	timestamp = JsonOutputStream(os.path.join(self.odir, 'timestamp.json'))
	timestamp.write(TIMESTAMP)
	timestamp.close()
        # All genomes file
        genomes = list(self.getGenomes())
        genomes.sort(lambda a,b: cmp(a['name'],b['name']))
        if (self.sample) :
          genomes = genomes[0:2]
        allGenomes = JsonOutputStream(os.path.join(self.odir, 'allGenomes.json'), indent=1)
        for g in genomes:
            self.log(g['name'])
	    gDir = os.path.join(self.odir, g['filename'])
	    self.ensureDir(gDir)
            allGenomes.write(g)
            # Init chromosome file
            chromosomes = JsonOutputStream(os.path.join(gDir, 'chromosomes.json'), indent=1)
            # Init feature file
            features = JsonOutputStream(os.path.join(gDir, 'features.json'), indent=1)
            # Process genome one chromosome at a time
            for ccount, c in enumerate(self.getChromosomes(g['name'])):
                if ccount > 2 :
                  continue
                # Write chromosome record
                self.log(c['name'])
                chromosomes.write(c)
                # Write features on this chromosome
                for f in self.getGenes(g['name'], c['name']):
                  # tp = 'pseudogene' if 'pseudo' in f['sotype'] else 'gene'
                  features.write(f)
                # Write models for this chromosome
                models = JsonOutputStream(os.path.join(gDir, 'models.chr%s.json'%c['name']), indent=1)
                for m in self.getModels(g['name'], c['name']):
                  models.write(m)
                models.close()
	    #
            chromosomes.close()
            features.close()
        # end for
        allGenomes.close()

class JsonOutputStream:
  def __init__ (self, ostream, indent=None, separators=(',', ': ')) :
    if type(ostream) == types.StringType:
      self.ostream = open(ostream, 'w')
    else:
      self.ostream = ostream
    self.count = 0
    self.separators = separators
    self.indent = indent
    self.ostream.write('[')

  def write (self, nxt) :
    if self.count:
      self.ostream.write(self.separators[0])
    self.ostream.write(json.dumps(nxt, indent=self.indent, separators=self.separators))
    self.count += 1

  def close (self):
    self.ostream.write(']')
    self.ostream.close()

if __name__ == "__main__":
    opts = getOpts()
    DataGetter(opts).main()

