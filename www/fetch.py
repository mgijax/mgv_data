
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
import time
import argparse
from urllib.request import urlopen
import json
import string
import cgi
import cgitb
import subprocess

# -----------------------------------------
TABIX=os.environ["TABIX"]
SAMTOOLS=os.environ["SAMTOOLS"]

# -----------------------------------------
DEFAULT_LINE_LEN = 60
DATA_DIR="."
MAX_DESCRIPTORS=4000
MAX_LENGTH=100000000

# -----------------------------------------
def chunkString (s, n) :
    return [s[i:i+n] for i in range(0, len(s), n)]

# -----------------------------------------
def defaultHeader (desc) :
    return ">default"

# -----------------------------------------
def getSequenceFromFaidx (desc) :
    arrayify = lambda x : x if type(x) is list else [x]
    gurl = desc["genomeUrl"]
    gdir = gurl.replace("/"," ").split()[-1]
    path = "%s/%s/assembly.fasta.gz" % (DATA_DIR, gdir)

    desc["chromosome"]
    starts = arrayify(desc["start"])
    lengths = arrayify(desc["length"])
    args = [ "%s:%d-%d" % (desc["chromosome"], s, s+lengths[i]-1) for i,s in enumerate(starts) ]
    #
    command = [SAMTOOLS, "faidx", path] + args
    r = subprocess.check_output(command)
    r = ''.join(filter(lambda l: not l.startswith('>'), r.decode('utf8').split('\n')))
    return r
  
# -----------------------------------------
def getSequenceFromAssembly(desc):
    seq = getSequenceFromFaidx(desc)
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
# Args:
#   seqId - string or list of strings - seqIds of desired sequences in curie format.
#   A curie is an id with a recognized prefix such as 'ENSEMBL' or 'RefSeq'. Example curies:
#       ENSEMBL:MGP_CASTEiJ_P0038053
#       RefSeq:FBtr0078736
#       RefSeq:NM_001122733
#       UniProt:Q9D404
def getSequenceFromSeqfetch (desc) :
    # Source names recognized by seqfetch
    #     swissprot trembl sptrembl genbank refseq ensembl_mus_cdna ensembl_mus_prot
    # Note that ensembl_mus_cdna and ensembl_mus_prot are really the same thing and actually implement
    # generic Ensembl sequence retrieval. Here we'll use ensembl_mus_cdna.
    #
    seqfetchBaseUrl = "http://www.informatics.jax.org/seqfetch/tofasta.cgi"
    # maps curie prefix to name to use in seqfetch requests
    PREFIXMAP = {
        "ensembl" : "ensembl_mus_cdna",
        "refseq" : "genbank",
        "uniprot" : "swissprot",
    }
    if "seqId" in desc:
        if type(desc["seqId"]) is str:
            seqIds = [desc["seqId"]]
        else:
            seqIds = desc["seqId"]
    else:
        raise RuntimeError("No seqId ids found in descriptor: " + str(desc))
    seqIds.sort()
    seqfetchDescrs = []
    for c in seqIds:
        if not ":" in c:
            raise RuntimeError("ID is not a curie: " + str(c))
        prefix, base = c.split(":", 1)
        np = PREFIXMAP[prefix.lower()]
        seqfetchDescrs.append("%s!%s!!!!!" % (np, base))
    seqfetchArgs = "&".join([ "seq%s=%s" % (i,d) for (i,d) in enumerate(seqfetchDescrs) ])
    url = seqfetchBaseUrl + '?' + seqfetchArgs
    with urlopen(url) as fd:
        seqs = fd.read()
        seqs = seqs.decode('utf-8')
    return seqs

# -----------------------------------------
def validateSequenceOptions (opts) :
    if not opts.descriptors:
        error("No descriptors.")
    ndescs = len(opts.descriptors)
    if ndescs > MAX_DESCRIPTORS:
        error("Too many descriptors: %d" % ndescs)
    tLength = 0
    for d in opts.descriptors:
        if "length" in d:
            if type(d["length"]) is list:
                tLength += sum(d["length"])
            else:
                tLength += d["length"]
    if tLength > MAX_LENGTH:
        error("Total length too big: %d" % tLength)

# -----------------------------------------
def doSequences (opts) :
    #
    validateSequenceOptions(opts)
    #
    print ('Content-Type: text/plain')
    if "filename" in opts:
        print(('Content-Disposition: attachment; filename = "%s"' % opts.filename))
    print ("")
    for d in opts.descriptors:
        if "seqId" in d:
            s = getSequenceFromSeqfetch(d)
        else:
            s = getSequenceFromAssembly(d)
        sys.stdout.write(s)
        sys.stdout.write('\n')

# -----------------------------------------
def getFeaturesFromTabix (desc) :
    gpath = desc["genome"]
    track = desc["track"]
    path = "%s/%s/%s.gff.gz" % (DATA_DIR, gpath, track)
    regions = desc["regions"]
    command = [TABIX, "--separate-regions", path] + regions.split()
    #
    r = subprocess.check_output(command)
    r = r.decode("utf8")
    h = '#genome=%s track=%s\n' % (gpath, track)
    return h+r

# -----------------------------------------
def doFeatures (opts) :
    if not opts.descriptors:
        error("No descriptors.")
    print ('Content-Type: text/plain')
    if "filename" in opts:
        print(('Content-Disposition: attachment; filename = "%s"' % opts.filename))
    print ("")
    for d in opts.descriptors:
        fs = getFeaturesFromTabix(d)
        sys.stdout.write(fs)

# -----------------------------------------
def doHomology (opts) :
    print ('Content-Type: text/plain')
    print ('')
    txid = opts.taxonid
    path = "%s/homologies/%s.tsv" % (DATA_DIR, txid)
    with open(path, 'r') as fd:
        sys.stdout.write(fd.read())

# -----------------------------------------
def doMetadata (opts) :
    metadata = []

    for name in os.listdir(DATA_DIR):
        gfile = os.path.join(DATA_DIR, name, "index.json")
        if not os.path.exists(gfile):
            continue
        with open(gfile, 'r') as fd:
            gcfg = json.load(fd)
        #
        metadata.append(gcfg)
        #
        # Add a timestamp, set to file's modification date.
        #
        filestats = os.stat(gfile)
        gcfg['timestamp'] = time.asctime(time.localtime(filestats.st_mtime))
        if gcfg['type'] == 'genome':
            #
            # Initialize list of chromosome names and lengths
            #
            chrFile = os.path.join(DATA_DIR, name, 'assembly.fasta.gz.fai')
            if os.path.exists(chrFile):
                with open(chrFile, 'r') as fd:
                    lines = fd.read().split('\n')[:-1]
                    rows = [ line.split('\t') for line in lines ]
                    chrs = [ { 'name':r[0], 'length':int(r[1]) } for r in rows ]
                    gcfg['chromosomes'] = chrs
            else:
                # Requires an indexed assembly to initialize the chromosomes assay. FIXME ?
                gcfg['chromosomes'] = []
        # end if
    # end for
    print('Content-Type: application/json')
    print('')
    print(json.dumps(metadata,indent=2))

# -----------------------------------------
def error (message) : 
   print ('Content-Type: text/plain')
   print ('')
   print ('ERROR: ' + message)
   sys.exit(1)

# -----------------------------------------
def getFormParameters (opts) :
    form = cgi.FieldStorage()
    if "datatype" in form:
        opts.datatype = form["datatype"].value
    if "descriptors" in form:
        opts.descriptors = form["descriptors"].value
    if "taxonid" in form:
        opts.taxonid = form["taxonid"].value
    if "filename" in form:
        opts.filename = form["filename"].value
    return opts

# -----------------------------------------
def getOptions () :
    parser = argparse.ArgumentParser(description="Get sequence slices from genome assembly sequences.")
    #
    parser.add_argument(
        "--datatype",
        dest="datatype",
        default="gff",
        choices=["fasta","gff","metadata","homology"],
        help="What to fetch. If fasta or gff, specific regions are specified in the descriptors argument. If metadata, returns object containing metadata for each known genome (descriptors, if provided, are ignored.)")

    parser.add_argument(
        "--descriptors",
        dest="descriptors",
        help="Descriptors, for datatype=fasta or gff.")

    parser.add_argument(
        "--taxonid",
        dest="taxonid",
        help="Taxon id, for datatype=homology.")

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
    #
    if opts.descriptors:
        if opts.descriptors == "TEST_FASTA":
            opts.datatype = "fasta"
            opts.descriptors = TEST_FASTA
        elif opts.descriptors == "TEST_GFF":
            opts.datatype = "gff"
            opts.descriptors = TEST_GFF
        else:
            print(opts.descriptors)
            opts.descriptors = json.loads(opts.descriptors)

    #
    return opts

# -----------------------------------------
def main () :
    opts = getOptions()
    global DATA_DIR
    DATA_DIR = opts.dataDirectory
    if opts.datatype == "metadata":
        doMetadata(opts)
    elif opts.datatype == "homology":
        doHomology(opts)
    elif opts.datatype == "fasta":
        doSequences(opts)
    elif opts.datatype == "gff":
        doFeatures(opts)
    else:
        raise RuntimeError("Unknown datatype value: " + opts.datatype)

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

TEST_GFF = [{
    "genome": "mus_musculus_aj",
    "track" : "models",
    "regions" : "1:123456790-124456790 2:223456790-224456790"
},{
    "genome": "mus_musculus_grcm39",
    "track" : "models.genes",
    "regions" : "18 19"
}]

TEST_FASTA = [{
    "genome": "mus_musculus_aj",
    "track": "assembly",
    "regions": "1:123456790-123456889",
    "header": ">test.0",
    #"lineLength" : 40,
    #"reverseComplement": False,
    #"translate": False,
}]

xTEST_FASTA = [{

"seqId" : [ "ENSEMBL:MGP_CASTEiJ_T0038053", "ENSEMBL:MGP_CASTEiJ_P0038053"]

}, {

"seqId" : "RefSeq:NM_001145293"

}, {

"header": ">test.0",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.0 line length = 40",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": False,
"translate": False,
"lineLength" : 40,

}, {

"header": ">test.0 rc",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": True,
"translate": False,

}, {

"header": ">test.0 xl",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": False,
"translate": True,

}, {

"header": ">test.0 rc xl",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 123456790,
"length": 100,
"reverseComplement": True,
"translate": True,

}, {

"header": ">test.1 pt1",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": [123456790],
"length": [100],
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 pt2",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 123458101,
"length": 102,
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 pt3",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": [123459222],
"length": [88],
"reverseComplement": False,
"translate": False,

}, {

"header": ">test.1 concatenated",
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": [123456790, 123458101, 123459222],
"length": [100, 102, 88],
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUSE00000702887 1:10038217-10038344 from File",
"reverseComplement": False,
"translate": False,
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": 10038217,
"length": 128,
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUST00000134716 from file",
"reverseComplement": False,
"translate": False,
"genome": "mus_musculus_aj",
"genomeUrl" : "mus_musculus_aj",
"chromosome": "1",
"start": [10039823,10045076,10047424,10058720,10059824],
"length": [57,109,103,104,168],
"reverseComplement": False,
"translate": False,

}, {

"header": ">ENSMUSE00001282101 1:10034008-10034136(-) from File",
"genome": "mus_musculus_aj",
"genomeUrl": "mus_musculus_aj",
"chromosome": "1",
"start": 10034008,
"length": 129,
"reverseComplement": True,
"translate": False,

}, {

"header": ">ENSMUST00000186528 from file",
"genome": "mus_musculus_aj",
"genomeUrl": "mus_musculus_aj",
"chromosome": "1",
"start": [10027198,10030599,10032339,10033258,10034008,10034929,10037662],
"length": [54,112,86,66,129,235,280],
"reverseComplement": True,
"translate": False,
"lineLength": 40

}, {

"header": ">ENSMUSP00000000834 from file",
"genome": "mus_musculus_aj",
"genomeUrl": "mus_musculus_aj",
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
amino_acids = [a.strip().split() for a in amino_acids_s.strip().split('\n')]
# create map from 3-letter AA code to single letter code
aaShort2Letter = dict([(r[1],r[2]) for r in amino_acids])
# create table of the genetic code from the string
genetic_code_t = [c.strip().split() for c in genetic_code_s.strip().split('\n')]
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
  residues = ''.join([aaShort2Letter.get(genetic_code.get(c, ''), '') for c in codons])
  return residues.split('X')[0] #truncate at first stop codon if there is one

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
