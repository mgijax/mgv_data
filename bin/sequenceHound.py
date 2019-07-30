#!/anaconda3/bin/python
#
# CGI for retrieving sequences in fasta format from multiple InterMine source and returniing
# the result in a single response. 
# Input is a list of descriptors, one per sequence.
# Response is a Fasta file.
#
# There are two ways to specify a sequence: as a slice of a chromosome, and as the sequence
# of a database object (like a CDS).
#
# -----------------------------------------
import sys
import os
import argparse
import urllib.request
import json
import string
import cgi
import cgitb
#cgitb.enable()

# -----------------------------------------
VALID_PREFIXES = [
  'http://',
  'https://'
]
VALID_WSROOTS = [
  'www.mousemine.org/mousemine/service/',
  'www.humanmine.org/humanmine/service/'
]
VALID_ENDPOINTS = [
  'query/results/fasta?',
  'sequence?'
]

# -----------------------------------------
def chunkString (s, n) :
  return [s[i:i+n] for i in range(0, len(s), n)]

# -----------------------------------------
def defaultHeader (desc) :
  desc['rc'] = 'reverse complement ' if desc['reverseComplement'] else ''
  desc['tr'] = 'translated ' if desc['translate'] else ''
  if 'chr' in desc:
    return '>%(genome)s::%(chr)s:%(start)d..%(end)d (%(rc)s%(tr)s%(type)s)' % desc
  else:
    return '>%(genome)s::%(ID)s (%(rc)s%(tr)s%(type)s)' % desc

# -----------------------------------------
def splitPrefix (s, possibles) :
  for p in possibles:
    if s.startswith(p):
      return (p, s[len(p):])
  return (None, s)

# -----------------------------------------
def validateUrl (url) :
  protocol, rest = splitPrefix(url, VALID_PREFIXES)
  if not protocol: return False
  #
  webservice, rest = splitPrefix(rest, VALID_WSROOTS)
  if not webservice: return False
  #
  endpoint, rest = splitPrefix(rest, VALID_ENDPOINTS)
  if not endpoint: return False
  #
  return True

# -----------------------------------------
# A mine sequence descriptor has these fields:
#   url (required) the fully formed url to a mine service endpoint
#   reverseComplement (optional) if truthy, reverse complements the returned sequence
#   translate (optional) if truthy, translates the sequence (must be a coding sequence)
#   header (optional) if provided, used as the header string for the sequence.
#        Otherise header is generated automatically.
def getMineSequence (desc) :
  if not validateUrl(desc['url']):
    return ''
  #
  fd = urllib.request.urlopen(desc['url'])
  data = fd.read().decode('utf-8')
  fd.close()
  #
  if data.startswith(">") :
    seq = ''.join(data.split('\n', 1)[1].split())
  elif data.startswith("{") :
    obj = json.loads(data)
    seq = obj['features'][0]['seq']
  #
  return seq

# -----------------------------------------
# Returns a sequence of characters from an open file
# Args:
#   fd - open file descriptor (as from os.open)
#   start - zero-based character position
#   length - number of characters to read
def slice (fd, start, length):
  os.lseek(fd, start, os.SEEK_SET) 
  return os.read(fd, length).decode('utf-8')

# -----------------------------------------
# 
def getFileSequence (desc) :
  # compose name of the file containing the chromosome sequence
  gname = desc["genome"]
  gfname = gname.replace("/","").lower()
  path = "./output/%s/sequences/%s" % (gfname, desc["chromosome"])
  fd = os.open(path, os.O_RDONLY)
  #
  starts = desc["start"]
  if type(starts) is int:
    starts = [starts]
  #
  lengths = desc["length"]
  if type(lengths) is int:
    lengths = [lengths]
  #
  seqs = []
  for i,s in enumerate(starts):
    l = lengths[i]
    seqs.append(slice(fd, s, l))
  seq = ''.join(seqs).lower()
  #
  return seq

# -----------------------------------------
def getSequence(desc):
  if "url" in desc:
     seq = getMineSequence(desc)
  else:
     seq = getFileSequence(desc)
  doRC = desc.get('reverseComplement', True)
  if desc.get('reverseComplement', False):
    seq = reverseComplement(seq)
  if desc.get('translate', False):
    seq = translate(seq)
  seq = '\n'.join(chunkString(seq, desc.get('lineLength', 60)))
  hdr = desc.get('header', None)
  if hdr is None:
      hdr = defaultHeader(desc)
  if not hdr.startswith('>'):
    hdr = '>' + hdr
  return hdr + '\n' + seq + '\n'

# -----------------------------------------
def doSequences (descs) :
  for d in descs:
    sys.stdout.write(getSequence(d))
    sys.stdout.write('\n')

# -----------------------------------------
def getParameters () :
  form = cgi.FieldStorage()
  params = {
    "descriptors" : None,
    "filename": "mgv.download.fa"
  }
  if "descriptors" in form:
      params["descriptors"] = json.loads(form["descriptors"].value)
  if "filename" in form:
      params["filename"] = form["filename"].value
  return params

# -----------------------------------------
def getCmdLineOptions () :
  parser = argparse.ArgumentParser(description="Get sequence slices from genome assembly sequences.")
  #
  parser.add_argument(
    "-T",
    dest="doTests",
    default=False,
    action="store_true",
    help="Run with test inputs.")
  clopts = parser.parse_args()
  params = {}
  if clopts.doTests:
    params["descriptors"] = TESTDATA
    params["filename"] = "testResults.fa"
  return params

# -----------------------------------------
def main () :
  doCGI = sys.argv[0].endswith('.cgi')
  if doCGI:
      params = getParameters()
  else:
      params = getCmdLineOptions()
  print ('Content-Type: text/x-fasta')
  print ('Content-Disposition: attachment; filename = "%s"' % params["filename"])
  print ()
  doSequences(params["descriptors"])

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

TESTDATA = [{

"header": ">test.0",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 123456789,
"length": 100,
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.0 line length = 40",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 123456789,
"length": 100,
"reverseComplement": False,
"translate": False,
"lineLength" : 40,

}, {

"header": ">test.0 rc",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 123456789,
"length": 100,
"reverseComplement": True,
"translate": False,

}, {

"header": ">test.0 xl",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 123456789,
"length": 100,
"reverseComplement": False,
"translate": True,

}, {

"header": ">test.0 rc xl",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 123456789,
"length": 100,
"reverseComplement": True,
"translate": True,

}, {

"header": ">test.1 pt1",
"genome": "C57BL/6J",
"chromosome": "1",
"start": [123456789],
"length": [100],
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 pt2",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 123458100,
"length": 102,
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 pt3",
"genome": "C57BL/6J",
"chromosome": "1",
"start": [123459221],
"length": [88],
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 concatenated",
"genome": "C57BL/6J",
"chromosome": "1",
"start": [123456789, 123458100, 123459221],
"length": [100, 102, 88],
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUSE00000702887 from MouseMine",
"reverseComplement": False,
"translate": False,
"url": "http://www.mousemine.org/mousemine/service/query/results/fasta?query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Exon.primaryIdentifier%22%3E%3Cconstraint%20path%3D%22Exon.primaryIdentifier%22%20op%3D%22ONE%20OF%22%3E%3Cvalue%3EENSMUSE00000702887%3C%2Fvalue%3E%3C%2Fconstraint%3E%3C%2Fquery%3E&view=Exon.gene.canonical.primaryIdentifier"

}, {

"header": ">ENSMUSE00000702887 1:10038217-10038344 from File",
"reverseComplement": False,
"translate": False,
"genome": "C57BL/6J",
"chromosome": "1",
"start": 10038216,
"length": 128,
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUST00000134716 from file",
"reverseComplement": False,
"translate": False,
"genome": "C57BL/6J",
"chromosome": "1",
"start": [10039822,10045075,10047423,10058719,10059823],
"length": [57,109,103,104,168],
"reverseComplement": False,
"translate": False,

}, { 

"header": ">ENSMUST00000134716 from MouseMine",
"reverseComplement": False,
"translate": False,
"url": "http://www.mousemine.org/mousemine/service/query/results/fasta?query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Transcript.primaryIdentifier%22%3E%3Cconstraint%20path%3D%22Transcript.primaryIdentifier%22%20op%3D%22ONE%20OF%22%3E%3Cvalue%3EENSMUST00000134716%3C%2Fvalue%3E%3C%2Fconstraint%3E%3C%2Fquery%3E&view=Transcript.gene.canonical.primaryIdentifier"

}, {

"header": ">ENSMUSE00001282101 1:10034008-10034136(-) from File",
"genome": "C57BL/6J",
"chromosome": "1",
"start": 10034007,
"length": 129,
"reverseComplement": True,
"translate": False,

}, {

"header": ">ENSMUSE00001282101 from MouseMine",
"reverseComplement": False,
"translate": False,
"url": "http://www.mousemine.org/mousemine/service/query/results/fasta?query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Exon.primaryIdentifier%22%3E%3Cconstraint%20path%3D%22Exon.primaryIdentifier%22%20op%3D%22ONE%20OF%22%3E%3Cvalue%3EENSMUSE00001282101%3C%2Fvalue%3E%3C%2Fconstraint%3E%3C%2Fquery%3E&view=Exon.gene.canonical.primaryIdentifier"

}, {

"header": ">ENSMUST00000186528 from file",
"genome": "C57BL/6J",
"chromosome": "1",
"start": [10027197,10030598,10032338,10033257,10034007,10034928,10037661],
"length": [54,112,86,66,129,235,280],
"reverseComplement": True,
"translate": False,
"lineLength": 60

}, {

"header": ">ENSMUST00000186528 from MouseMine",
"reverseComplement": False,
"translate": False,
"url": "http://www.mousemine.org/mousemine/service/query/results/fasta?query=%3Cquery%20model%3D%22genomic%22%20view%3D%22Transcript.primaryIdentifier%22%3E%3Cconstraint%20path%3D%22Transcript.primaryIdentifier%22%20op%3D%22ONE%20OF%22%3E%3Cvalue%3EENSMUST00000186528%3C%2Fvalue%3E%3C%2Fconstraint%3E%3C%2Fquery%3E&view=Transcript.gene.canonical.primaryIdentifier"

}, {

"header": ">ENSMUSP00000000834 from file",
"genome": "C57BL/6J",
"chromosome": "1",
"start": [161781575,161782947,161787104,161787943],
"length": [395,57,46,342],
"reverseComplement": True,
"translate": True,

}, {

"header": ">ENSMUSP00000000834 from MouseMine",
"reverseComplement": False,
"translate": True,
"url": "http://www.mousemine.org/mousemine/service/query/results/fasta?query=%3Cquery%20model%3D%22genomic%22%20view%3D%22CDS.primaryIdentifier%22%3E%3Cconstraint%20path%3D%22CDS.primaryIdentifier%22%20op%3D%22ONE%20OF%22%3E%3Cvalue%3EENSMUSP00000000834%3C%2Fvalue%3E%3C%2Fconstraint%3E%3C%2Fquery%3E&view=CDS.gene.canonical.primaryIdentifier"

}]

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

# maps a base to its complement
base_complement = { 
  'a' : 't',
  't' : 'a',
  'c' : 'g',
  'g' : 'c',
  'n' : 'n',
  
  'A' : 'T',
  'T' : 'A',
  'C' : 'G',
  'G' : 'C',
  'N' : 'N' 
}

# table of amino acids.
amino_acids_s = ''' 
  Alanine Ala     A   
  Arginine        Arg     R   
  Asparagine      Asn     N   
  Aspartate       Asp     D   
  Cysteine        Cys     C   
  Glutamate       Glu     E   
  Glutamine       Gln     Q   
  Glycine Gly     G   
  Histidine       His     H   
  Isoleucine      Ile     I   
  Leucine Leu     L   
  Lysine  Lys     K   
  Methionine      Met     M   
  Phenylalanine   Phe     F   
  Proline Pro     P   
  Selenocysteine  Sec     U   
  Serine  Ser     S   
  Threonine       Thr     T   
  Tryptophan      Trp     W   
  Tyrosine        Tyr     Y   
  Valine  Val     V   
  Stop      Stop    X   
'''

# The genetic code. Map each codon to the abbrev of its translated residue.
genetic_code_s = '''
  UUU     Phe
  UUC     Phe
  UUA     Leu
  UUG     Leu
  UCU     Ser
  UCC     Ser
  UCA     Ser
  UCG     Ser
  UAU     Tyr
  UAC     Tyr
  UAA     Stop
  UAG     Stop
  UGU     Cys
  UGC     Cys
  UGA     Stop
  UGG     Trp
  CUU     Leu
  CUC     Leu
  CUA     Leu
  CUG     Leu
  CCU     Pro
  CCC     Pro
  CCA     Pro
  CCG     Pro
  CAU     His
  CAC     His
  CAA     Gln
  CAG     Gln
  CGU     Arg
  CGC     Arg
  CGA     Arg
  CGG     Arg
  AUU     Ile
  AUC     Ile
  AUA     Ile
  AUG     Met
  ACU     Thr
  ACC     Thr
  ACA     Thr
  ACG     Thr
  AAU     Asn
  AAC     Asn
  AAA     Lys
  AAG     Lys
  AGU     Ser
  AGC     Ser
  AGA     Arg
  AGG     Arg
  GUU     Val
  GUC     Val
  GUA     Val
  GUG     Val
  GCU     Ala
  GCC     Ala
  GCA     Ala
  GCG     Ala
  GAU     Asp
  GAC     Asp
  GAA     Glu
  GAG     Glu
  GGU     Gly
  GGC     Gly
  GGA     Gly
  GGG     Gly
'''

# create table of amino acids from the string
amino_acids = map(lambda a: a.strip().split(), amino_acids_s.strip().split('\n'))
# create map from 3-letter AA code to single letter code
aaShort2Letter = dict([(r[1],r[2]) for r in amino_acids])
# create table of the genetic code from the string
genetic_code_t = map(lambda c: c.strip().split(), genetic_code_s.strip().split('\n'))
# create map from codon to residue
genetic_code = dict([ (r[0],r[1]) for r in genetic_code_t ])
# create character mapping table for doing complements
compTable = ''.maketrans('actgACTG', 'tgacTGAC')

def translate (cds) :
  cds = cds.upper().replace('T', 'U')
  codons = [ cds[i:i+3] for i in range(0, len(cds),3) ]
  residues = ''.join(map(lambda c: aaShort2Letter.get(genetic_code.get(c, ''), ''), codons))
  return residues

def complement (dna) :
  return ''.join([base_complement.get(b, b) for b in dna])

def reverseComplement (dna) :
  return complement(dna)[::-1]

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

# GO!!
main()
