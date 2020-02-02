
# fetch.py
#
# CGI for retrieving an arbitrary set of sequences from the available genomes.
#
# Input is a list of descriptors, one per sequence.
# Each descriptor specifies a genome, a chromosome, and one or more start/length pairs.
#
# When multiple start/length pairs are specified, the pieces are concatenated.
# A descriptor may specify that the sequence should be reverse complemented.
# A descriptor may specify that the sequence should be translated.
# (The order of things is: get the pieces, concatenate, reverse complement, translate.)
# A descriptor specifies the header to be returned with the sequence.
#
# Response is a Fasta file.
#
# -----------------------------------------
import sys
import os
import argparse
try:
  from urllib.request import urlopen
except:
  from urllib import urlopen
import json
import string
import cgi
import cgitb
#cgitb.enable()

# -----------------------------------------
DEFAULT_LINE_LEN = 60
DATA_DIR="."

# -----------------------------------------
def chunkString (s, n) :
  return [s[i:i+n] for i in range(0, len(s), n)]

# -----------------------------------------
def defaultHeader (desc) :
  return ">default"

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
def getSequenceFromFile (desc) :
  # compose name of the file containing the chromosome sequence
  # The directory is the last component of the url:
  gurl = desc["genomeUrl"]
  gdir = gurl.replace("/"," ").split()[-1]
  path = "%s/%s/sequences/%s.txt" % (DATA_DIR, gdir, desc["chromosome"])
  # 
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
    seqs.append(slice(fd, s-1, l))
  seq = ''.join(seqs).lower()
  #
  return seq

# -----------------------------------------
def getSequence(desc):
  seq = getSequenceFromFile(desc)
  if desc.get('reverseComplement', False):
    seq = reverseComplement(seq)
  if desc.get('translate', False):
    seq = translate(seq)
  seq = '\n'.join(chunkString(seq, desc.get('lineLength', DEFAULT_LINE_LEN)))
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
def getFormParameters (opts) :
  form = cgi.FieldStorage()
  if "descriptors" in form:
      opts.descriptors = json.loads(form["descriptors"].value)
  if "filename" in form:
      opts.filename = form["filename"].value
  return opts

# -----------------------------------------
def getOptions () :
  parser = argparse.ArgumentParser(description="Get sequence slices from genome assembly sequences.")
  #
  parser.add_argument(
    "--descriptors",
    dest="descriptors",
    action="append",
    help="Descriptors.")

  parser.add_argument(
    "--dir",
    dest="dataDirectory",
    default=".",
    help="Path to data directory.")

  parser.add_argument(
    "--cgi",
    dest="doCGI",
    action="store_true",
    default=False,
    help="Run as a CGI script.")

  opts = parser.parse_args()
  if opts.doCGI:
    opts = getFormParameters(opts)
  if not opts.descriptors:
    opts.descriptors = TESTDATA
  return opts

# -----------------------------------------
def main () :
  opts = getOptions()
  global DATA_DIR
  DATA_DIR = opts.dataDirectory
  print ('Content-Type: text/plain')
  #print ('Content-Type: text/x-fasta')
  if "filename" in opts:
     print ('Content-Disposition: attachment; filename = "%s"' % opts.filename)
  print ("")
  doSequences(opts.descriptors)

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

TESTDATA = [{

"header": ">test.0",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.0 line length = 40",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": False,
"translate": False,
"lineLength" : 40,

}, {

"header": ">test.0 rc",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": True,
"translate": False,

}, {

"header": ">test.0 xl",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": False,
"translate": True,

}, {

"header": ">test.0 rc xl",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": True,
"translate": True,

}, {

"header": ">test.1 pt1",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": [123456790],
"length": [100],
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 pt2",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 123458101,
"length": 102,
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 pt3",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": [123459222],
"length": [88],
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 concatenated",
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": [123456790, 123458101, 123459222],
"length": [100, 102, 88],
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUSE00000702887 1:10038217-10038344 from File",
"reverseComplement": False,
"translate": False,
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": 10038217,
"length": 128,
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUST00000134716 from file",
"reverseComplement": False,
"translate": False,
"genome": "mus_musculus",
"genomeUrl" : "mus_musculus",
"chromosome": "1",
"start": [10039823,10045076,10047424,10058720,10059824],
"length": [57,109,103,104,168],
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUSE00001282101 1:10034008-10034136(-) from File",
"genome": "mus_musculus",
"genomeUrl": "mus_musculus",
"chromosome": "1",
"start": 10034008,
"length": 129,
"reverseComplement": True,
"translate": False,

}, {

"header": ">ENSMUST00000186528 from file",
"genome": "mus_musculus",
"genomeUrl": "mus_musculus",
"chromosome": "1",
"start": [10027198,10030599,10032339,10033258,10034008,10034929,10037662],
"length": [54,112,86,66,129,235,280],
"reverseComplement": True,
"translate": False,
"lineLength": 40

}, {

"header": ">ENSMUSP00000000834 from file",
"genome": "mus_musculus",
"genomeUrl": "mus_musculus",
"chromosome": "1",
"start": [161781576,161782948,161787105,161787944],
"length": [395,57,46,342],
"reverseComplement": True,
"translate": True,

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
try:
  compTable = ''.maketrans('actgACTG', 'tgacTGAC')
except:
  import string
  compTable = string.maketrans('actgACTG', 'tgacTGAC')


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
if __name__ == '__main__':
  main()
