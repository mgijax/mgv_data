
# fetch.py
#
# CGI for retrieving an features, sequences, and metadata for a set of deployed genomes.
# For retrievig deployed homology data.
#
# FETCHING METADATA
#   Args:
#     datatype : metadata
#   Response:
#     Response is a JSON object containing metadata for all deployed genomes.
#     MIME type: application/json
#
# FETCHING FEATURES
#   Args:
#     datatype : gff
#     track : <trackname> E.g., models or models.genes
#     descriptors : <JSON encoded list> Each list element is an object with two fields:
#               genome : <genomepath> Path name for the genome e.g., mus_musculus_aj
#               regions : <regions> Argument to pass to tabix, e.g. "1:123456789-123498765"
#                       Specify multiple regions as space separated list. Specify whole chromosomes 
#                       by just naming the chromosome.
#     filename : <filename> to add as a Content-Disposition response header
#
#  Response:
#     Response is a GFF file of features 
#     MIME type: application/x-gff
#               
# FETCHING SEQUENCES
#  Args:
#     datatype : fasta
#     track : <trackname> E.g., assembly
#     descriptors : <JSON encoded list> Each list element describes one sequence to return.
#         Each descriptor has these fields:
#             genome : <genomepath> Path name for the genome, e.g., mus_musculus_aj
#             regions : <regions> Argument for faidx. One or more space separated regions, 
#                       each of the form <chr>:<start>-<end>. If multiple regions are specified, the
#                       sequences are concatenated to a single result sequence.
#             reverse : <boolean> If true, reverse complement the sequence. Default = false
#             translate : <boolean> If true, translate the sequence. Default = false
#             header : <string> If provided, used as the header line for the returned sequence.
#                       Otherwise a default header is generated.
#     filename : <filename> to add as a Content-Disposition response header
#
#  Response:
#     Response is a Fasta file of sequences, one per descriptor. 
#     There are hard limits on the number of descriptors and total size of the request.
#     MIME type: application/x-fasta
#
# FETCHING HOMOLOGY DATA
#  Args:
#    datatype : homology
#    taxonid : <taxonid> (the numeric part of an NBCI Taxon identifier, e.g. 9606 for human.
#  Response:
#    Returns orthology data for the given taxon id. Returns a table in TSV format with four columns:
#         geneid1, taxon1, geneid2, taxon2
#    In all rows, taxon1 is equal to the specified taxon
#    MIME type: text/tab-separated-values
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
import re

# -----------------------------------------
TABIX=os.environ["TABIX"]
SAMTOOLS=os.environ["SAMTOOLS"]
DEVMODE=os.environ["DEVMODE"]

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
    hdr = ">%s::%s" % (desc["genome"], "+".join(desc["regions"].split()))
    if desc.get("reverse", False):
        hdr += ",RC"
    if desc.get("translate", False):
        hdr += ",M"
    return hdr

# -----------------------------------------
def getSequenceFromFaidx (desc) :
    arrayify = lambda x : x if type(x) is list else [x]
    gpath = desc["genome"]
    track = desc["track"]
    path = "%s/%s/%s.fasta.gz" % (DATA_DIR, gpath, track)
    regions = desc["regions"]
    #
    command = [SAMTOOLS, "faidx", path] + regions.split()
    r = subprocess.check_output(command)
    r = ''.join(filter(lambda l: not l.startswith('>'), r.decode('utf8').split('\n')))
    return r
  
# -----------------------------------------
def getSequenceFromAssembly(desc):
    seq = getSequenceFromFaidx(desc)
    if desc.get('reverse', False):
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
def validateSequenceOptions (opts) :
    if not opts.descriptors:
        error("No descriptors.")
    ndescs = 0
    tLength = 0
    region_re = re.compile("^([^:]+):(\d+)-(\d+)$")
    for d in opts.descriptors:
        rs = d["regions"].split()
        ndescs += len(rs)
        if ndescs > MAX_DESCRIPTORS:
            error("Too many descriptors")
        for r in rs:
            m = region_re.match(r)
            if not m:
                error("Bad region descriptor.")
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            tLength += end - start + 1
            if tLength > MAX_LENGTH:
                error("Total request size too big.")

# -----------------------------------------
def doSequences (opts) :
    #
    validateSequenceOptions(opts)
    #
    print ('Content-Type: application/x-fasta')
    printCORS()
    if "filename" in opts:
        print(('Content-Disposition: attachment; filename = "%s"' % opts.filename))
    print ("")
    for d in opts.descriptors:
        if "seqId" in d:
            raise RuntimeError("Not implemented.")
        else:
            s = getSequenceFromAssembly(d)
        sys.stdout.write(s)

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
    print ('Content-Type: application/x-gff')
    printCORS()
    if "filename" in opts:
        print(('Content-Disposition: attachment; filename = "%s"' % opts.filename))
    print ("")
    for d in opts.descriptors:
        fs = getFeaturesFromTabix(d)
        sys.stdout.write(fs)

# -----------------------------------------
def doHomology (opts) :
    print ('Content-Type: text/tab-separated-values')
    printCORS()
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
    printCORS()
    print('')
    print(json.dumps(metadata,indent=2))

# -----------------------------------------
def error (message) : 
   print ('Content-Type: text/plain')
   printCORS()
   print ('')
   print ('ERROR: ' + message)
   sys.exit(1)

# -----------------------------------------
def printCORS () :
    if DEVMODE == "true":
        print('Access-Control-Allow-Origin: *')

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
        choices=["fasta","gff","metadata","homology"],
        help="What to fetch. If fasta or gff, specific regions are specified in the descriptors argument. If metadata, returns object containing metadata for each known genome (descriptors, if provided, are ignored.)")

    parser.add_argument(
        "--descriptors",
        dest="descriptors",
        help="Descriptors, for datatype=fasta or gff. Can also specify one of the test sets of descriptor by name. Choose one of: %s" % ', '.join(TESTS.keys())  )

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
        if opts.descriptors.startswith("TEST_"):
            tname = opts.descriptors
            if tname not in TESTS:
                error("Bad test name.")
            opts.datatype = TESTS[tname][0]
            opts.descriptors = TESTS[tname][1]
        else:
            try:
                opts.descriptors = json.loads(opts.descriptors)
            except:
                error("Bad descriptors.")
    #
    if not opts.datatype:
        parser.error("Please specify a datatype.")
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

TESTS = {
    "TEST_GFF" : ("gff",[{
        "genome": "mus_musculus_aj",
        "track" : "models",
        "regions" : "1:123456790-124456790 2:223456790-224456790"
    },{
        "genome": "mus_musculus_grcm39",
        "track" : "models.genes",
        "regions" : "18 19"
    }]),
    "TEST_FASTA_TOO_BIG" : ("fasta", [{
        "header": "Too big.",
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:1-123456889",
    }]),
    "TEST_FASTA" : ("fasta", [{
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:123456790-123456889",
    }, {
        "header": ">retrieve 100 bp, specify header",
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:133456790-133456889",
    }, {
        "header": ">retrieve 200 bp, joined",
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:123456790-123456889 1:133456790-133456889",
    }, {
        "header": ">retrieve 200 bp, set line length",
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:123456790-123456889 1:133456790-133456889",
        "lineLength" : 40,
    }, {
        "header": ">retrieve 200 bp, reverse complement",
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:123456790-123456889 1:133456790-133456889",
        "reverse": True,
        "translate": False,
    }, {
        "header": ">retrieve 200 bp, translate",
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:123456790-123456889 1:133456790-133456889",
        "reverse": False,
        "translate": True,
    }, {
        "genome": "mus_musculus_aj",
        "track": "assembly",
        "regions": "1:123456790-123456889 1:133456790-133456889",
        "reverse": True,
        "translate": True,

    #}, {
    #    "seqId" : [ "ENSEMBL:MGP_CASTEiJ_T0038053", "ENSEMBL:MGP_CASTEiJ_P0038053"]
    #}, {
    #    "seqId" : "RefSeq:NM_001145293"

    }])
}

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
