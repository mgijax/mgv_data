#
# downloadEnsemblGenome.py
#
# Downloads a genome annotation (GFF3) file or genome assembly (FASTA) file 
# from Ensembl for the specified organism and release version.
# Writes the UNcompressed data to stdout.
#
# Example:
#
# Download release 96 genome annotations for A/J and write to a GFF3 file:
#   python downloadEnsemblGenome.py -n mus_musculus_aj -r 96 -g models > aj.gff3
#
# Download the current release genome assembly for cow and pipe to another program.
#   python downloadEnsemblGenome.py -n bos_taurus -g assembly | python split.py ...
#
import sys
import os
import re
import argparse
import urllib.request
import gzip
from importGff3 import Gff3Importer, getArgs as gff3Args
import gff3lite

ENSEMBL_BASE = "ftp://ftp.ensembl.org/pub/"
GFF3_EXCLUDE = ["chromosome", "biological_region"]

############################################
## Opts
############################################

def getOpts () :
    parser = argparse.ArgumentParser()
    #
    parser.add_argument(
        "-n",
        dest="organism",
        required=True,
        help="Organism name used in ENSEMBL paths (eg, mus_musculus_aj)")
    #
    CHOICES = ['models', 'assembly']
    parser.add_argument(
        "-g",
        dest="get",
        required=True,
        choices=CHOICES,
        help="What to get. Required.")
    #
    parser.add_argument(
        "-x",
        dest="exclude",
        action="append",
        help="SO term of features to exclude from output. Repeatible. Default = " + str(GFF3_EXCLUDE))
    #
    parser.add_argument(
        "-r",
        dest="release",
        default="current",
        help="ENSEMBL release number (eg, 96). By default, gets the current release.")
    #
    parser.add_argument(
        "-c",
        dest="chromosomes",
        default='',
        help="List of chromosomes to include (comma separated, no spaces). By default, includes all chromosomes and contigs.")
    #
    parser.add_argument(
        "-C",
        dest="chrRegex",
        default='',
        help="Regular expression that chromosome names must match. By default, no matching is done.")
    #
    opts = parser.parse_args()
    if not opts.exclude:
      opts.exclude = GFF3_EXCLUDE
    if opts.chromosomes:
      opts.chromosomes = opts.chromosomes.split(',')
    opts.cset = set(opts.chromosomes)
    if opts.chrRegex:
      opts.chrRegex = re.compile(opts.chrRegex)
    return opts

############################################
## Utils
############################################

def getDirectoryList(url):
    fd = urllib.request.urlopen(url)
    dat = fd.read().decode('utf-8')
    fd.close()
    fnames = list(filter(lambda s: s.endswith('.gz'), dat.split()))
    return fnames

def ensureDir (d) :
    if not os.path.isdir(d):
        os.makedirs(d)
    return d

############################################
## FASTA
############################################

# Returns the URL for getting genome assembly (FASTA) file from Ensembl 
# for the specified organism and release.
#
def getFastaPath (organism, release):
    if release == "current":
        dirUrl = "%s/current_fasta/%s/dna" % (ENSEMBL_BASE,organism)
    else:
        dirUrl = "%s/release-%s/fasta/%s/dna" % (ENSEMBL_BASE,release,organism)
    fnames = list(filter(lambda n: n.endswith('.dna.toplevel.fa.gz'), getDirectoryList(dirUrl)))
    return (dirUrl + "/" + fnames[0], fnames[0])

# Downloads the genome assembly specified in opts, decompresses it, and writes to stdout.
#
def doFasta (opts) :
    (gpath, fname) = getFastaPath (opts.organism, opts.release)
    fd = urllib.request.urlopen(gpath)
    gfd = map(lambda b: b.decode('utf-8'), gzip.GzipFile(fileobj=fd))
    printIt = True
    for line in gfd:
      if line.startswith('>'):
        chrom = line[1:].split()[0]
        if opts.chromosomes:
          printIt = chrom in opts.chromosomes
        elif opts.chrRegex:
          printIt = opts.chrRegex.match(chrom) is not None
      #
      if printIt:
        print(line, end='')

############################################
## GFF3
############################################

# Returns the URL for getting genome annotation (GFF3) file from Ensembl for the specified
# organism and release.
#
def getGff3Path (organism, release="current"):
    if release == "current":
        dirUrl = "%s/current_gff3/%s" % (ENSEMBL_BASE,organism)
        suffix = '.gff3.gz'
        fnames = list(filter(lambda n: n.endswith(suffix) and len(n.split('.')) == 5, getDirectoryList(dirUrl)))
    else:
        dirUrl = "%s/release-%s/gff3/%s" % (ENSEMBL_BASE,release,organism)
        suffix = '.%s.gff3.gz' % release
        fnames = list(filter(lambda n: n.endswith(suffix), getDirectoryList(dirUrl)))
    ##
    return (dirUrl + "/" + fnames[0], fnames[0])

# Processes one line from the gff3, and either writes it to stdout or suppresses it.
# Filters by specified types (column 3) values and chromosomes.
def processGff3Line (line, opts):
    global lastLine
    if line.startswith('#'):
        if line.startswith('##sequence-region'):
            (tag, chrom, cstart, clength) = line[:-1].split()
            if opts.chromosomes and chrom not in opts.cset:
              return
            if opts.chrRegex and not opts.chrRegex.match(chrom):
              return
        elif line[:3] == "###":
            if line == lastLine:
              return
        print(line, end='')
    else:
        tokens = gff3lite.parseLine(line)
        if tokens[2] in opts.exclude:
          return
        if opts.chromosomes and tokens[0] not in opts.cset:
          return
        if opts.chrRegex and not opts.chrRegex.match(tokens[0]):
          return
        print(gff3lite.formatLine(tokens), end='')
    lastLine = line

   
# Downloads the genome annotation specified in the opts, decompresses it, and writes to stdout
#
def doGff3 (opts):
    (gpath, fname) = getGff3Path(opts.organism, opts.release)
    fd = urllib.request.urlopen(gpath)
    gfd = map(lambda b: b.decode('utf-8'), gzip.GzipFile(fileobj=fd))
    for line in gfd:
        processGff3Line(line, opts)

############################################
## Main
############################################

def main():
    opts = getOpts()
    if 'models' == opts.get:
      doGff3(opts)
    elif 'assembly' == opts.get:
      doFasta(opts)

#
main()
