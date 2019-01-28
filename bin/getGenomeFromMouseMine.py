#
# getGenomeFromMouseMine.py
#
# Generates a GFF3 file of gene models for the specified genomes.
# Writes to stdout.
#
import types
import sys
import urllib
import os.path
import json
import time
import itertools
import argparse
import gff3lite

# Web service address root for mousemine
MOUSEMINE="http://www.mousemine.org/mousemine/service"
#
TIMESTAMP=time.asctime(time.localtime(time.time()))

def getOpts () :
  parser = argparse.ArgumentParser(description='Dumps files for MGV from MouseMine.')
  parser.add_argument('-g', '--genome-name',
    dest='genomes', 
    action='append',
    default=[],
    help='Name of the genome(s) to dump. Repeatable. (default=dump all available genomes)')
  parser.add_argument('-d', '--output-dir',
    dest='outputDir', 
    default='./output',
    help='Output directory. Writes the gff3 file to this directory. Specify "-" to write to standard out. (default=%(default)s)')
  parser.add_argument('-t', '--timestamp',
    dest='timestamp', 
    default=TIMESTAMP,
    help='Timestamp to add to the file output file header. (default=current time).')
  parser.add_argument('--sample', dest='sample', 
    default=False,
    action='store_true',
    help='Generate sample output (default=%(default)s)')
  parser.add_argument(
    '-u', '--url',
    dest='url', 
    default=MOUSEMINE,
    help='Base web service URL. (default=%(default)s)')
  return parser.parse_args()

def doMouseMineQuery(wsurl, q):
    fmt = 'tab'
    url = '%s/query/results?format=%s&query=%s' % (wsurl,fmt,urllib.quote_plus(q))
    fd = urllib.urlopen(url)
    for line in fd:
        toks = line[:-1].split('\t')
        yield toks
    fd.close()

class GenomeDumper :
    def __init__ (self, opts) :
        self.gname = opts.gname
        self.outputDir = opts.outputDir
        self.url = opts.url
        self.sample = opts.sample
        self.lFd = sys.stderr
        self.ofd = sys.stdout
        self.timestamp = opts.timestamp
        self.taxonid = opts.taxonid

    def ensureDir(self, d):
        if not os.path.exists(d):
          os.makedirs(d)

    def open (self):
        if self.outputDir == '-':
          self.ofn = '<STDOUT>'
          self.ofd = sys.stdout
        else:
          self.ofn = os.path.join(self.outputDir, self.gname.replace('/','').lower() + '.gff3')
          self.ensureDir(self.outputDir)
          self.ofd = open(self.ofn, 'w')

    def log(self, s, NL='\n'):
        self.lFd.write(s)
        self.lFd.write(NL)

    def convertGenomeName (self, s):
        return s.lower().replace('/','')

    def convertStrand (self, s) :
        return '+' if s == 1 else '-'

    # Returns iterator over chromosomes for the given genome
    # Args:
    #    g (string) the name of the genome, eg, "A/J"
    # Returns:
    #    List of rows, one per chromosome of the given genome.
    #    Each row has a chromosome id and total chromosome length.
    #    The rows are sorted by chromosome id.
    def getChromosomes (self) :
        q = '''<query 
        model="genomic"
        view="Chromosome.primaryIdentifier Chromosome.length"
        sortOrder="Chromosome.symbol ASC"
        >
        <constraint path="Chromosome.strain.name" op="=" value="%s"/>
        </query>''' % self.gname

        for r in doMouseMineQuery(self.url, q):
            yield {
              'name': r[0],
              'length': int(r[1])
            }

    # Returns an iterator over genes on the specified chromosome for the specified genome.
    # Args:
    #    c (string) the chromosome, eg, "13"
    # Returns:
    #    Iterator that yields one row per gene in the specified genome.
    #    Each row has: strain name, gene id, chromosome location (chr, start, end, strand)
    def getGenes (self, c) :
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
        </query>''' % (self.gname, c)

        for r in doMouseMineQuery(self.url, q):
            attrs = {
              'ID': r[1],
              'transcripts': {}
            }
            if r[7] != '""':
                attrs['cID'] = r[7]
                attrs['symbol'] = r[8]
            gf = [
              r[3],
              'MGI',
              r[2],
              int(r[4]),
              int(r[5]),
              '.',
              self.convertStrand(r[6]),
              '.',
              attrs
            ]
            yield gf

    def getTranscripts (self, c):
      q = '''<query
      model="genomic"
      view="
      Gene.primaryIdentifier
      Gene.transcripts.primaryIdentifier
      Gene.transcripts.sequenceOntologyTerm.name
      Gene.transcripts.chromosomeLocation.start
      Gene.transcripts.chromosomeLocation.end
      Gene.transcripts.chromosomeLocation.strand
      "
      >
      <constraint path="Gene.strain.name" op="=" value="%s"/>
      <constraint path="Gene.chromosome.primaryIdentifier" op="=" value="%s"/>
      </query>
      ''' % (self.gname, c)
      for r in doMouseMineQuery(self.url, q):
          yield [
              c,
              'MGI',
              r[2],
              r[3],
              r[4],
              '.',
              self.convertStrand(r[5]),
              '.',
              { 'ID': r[1],
                'Parent': r[0],
                'exons' : []
              }
          ]

    def getExons (self, c):
      q = '''<query
      model="genomic"
      view="
      Gene.primaryIdentifier
      Gene.transcripts.primaryIdentifier
      Gene.transcripts.exons.primaryIdentifier
      Gene.transcripts.exons.chromosomeLocation.start
      Gene.transcripts.exons.chromosomeLocation.end
      Gene.transcripts.exons.chromosomeLocation.strand
      "
      >
      <constraint path="Gene.strain.name" op="=" value="%s"/>
      <constraint path="Gene.chromosome.primaryIdentifier" op="=" value="%s"/>
      </query>
      ''' % (self.gname, c)
      for r in doMouseMineQuery(self.url, q):
          yield [
              c,
              'MGI',
              'exon',
              r[3],
              r[4],
              '.',
              self.convertStrand(r[5]),
              '.',
              { 'ID': r[2],
                'Parent': r[1],
                'gene_id' : r[0]
              }
          ]

    def ensureDir(self, d):
        if not os.path.exists(d):
            os.makedirs(d)

    def getGenomeInfo (self):
        chrs = list(self.getChromosomes())
        if self.sample:
            chrs = chrs[0:2]
        self.genomeInfo = {
          "name" : self.gname,
          "taxonid" : self.taxonid,
          "chromosomes" : chrs,
          "timestamp" : TIMESTAMP
        }
        self.log('\nGenome: %s' % self.genomeInfo['name'])

    def doChromosome(self, c):
        # Assemble the models
        self.log('Chromosome %s ...' % c['name'], '')
        gid2g = {}
        for ig, g in enumerate(self.getGenes(c['name'])):
          attrs = g[8]
          gid = attrs['ID']
          gid2g[gid] = g
        self.log('%d genes ...' % (ig+1), '')
        for it, t in enumerate(self.getTranscripts(c['name'])):
          attrs = t[8]
          tid = attrs['ID']
          gid = attrs['Parent']
          gts = gid2g[gid][8]['transcripts']
          gts[tid] = t
        self.log('%d transcripts ...' % (it+1), '')
        for ie, e in enumerate(self.getExons(c['name'])):
          attrs = e[8]
          eid = attrs['ID']
          tid = attrs['Parent']
          gid = attrs['gene_id']
          ts = gid2g[gid][8]['transcripts']
          ts[tid][8]['exons'].append(e)
        self.log('%d exons' % (ie+1))
        # Output
        genes = gid2g.values()
        genes.sort(lambda a,b: a[3] - b[3])
        for g in genes:
          transcripts = g[8].pop('transcripts').values()
          self.ofd.write(gff3lite.formatLine(g))
          for t in transcripts:
            exons = t[8].pop('exons')
            self.ofd.write(gff3lite.formatLine(t))
            for e in exons:
              self.ofd.write(gff3lite.formatLine(e))
          # end for t
          self.ofd.write('###\n')
        # end for g

    def writeHeader(self):
        self.open()
        self.ofd.write(gff3lite.GFF3HEADER)
        self.ofd.write('##genome-name %s\n' % self.genomeInfo['name'])
        self.ofd.write('##taxonid %s\n' % self.genomeInfo['taxonid'])
        self.ofd.write('##date %s\n' % self.timestamp)
        for c in self.genomeInfo['chromosomes']:
          self.ofd.write('##sequence-region %s 1 %d\n' % (c['name'], c['length']))

    def main (self):
        self.getGenomeInfo()
        self.writeHeader()
        for c in self.genomeInfo['chromosomes']:
          self.doChromosome(c)

def getAvailableGenomes (url) :
    # returns strain name, taxonid pairs
    q='''
    <query model="genomic" view="Strain.name Strain.organism.taxonId" sortOrder="Strain.name asc">
      <constraint path="Strain" op="IN" value="Annotated strains"/>
      </query>
    '''
    return doMouseMineQuery(url, q)

if __name__ == "__main__":
    opts = getOpts()
    allGenomes = getAvailableGenomes(opts.url)
    if len(opts.genomes) == 0:
        genomes = list(allGenomes)
        if opts.sample:
            genomes = opts.genomes[0:2]
    else:
        genomes = []
        for gn in opts.genomes:
            gngs = filter(lambda g: gn == g[0], allGenomes)
            if len(gngs) == 0:
                raise RuntimeError('Genome %s not found.' % gn)
            genomes += gngs
    for g in genomes:
        opts.gname = g[0]
        opts.taxonid = g[1]
        GenomeDumper(opts).main()

