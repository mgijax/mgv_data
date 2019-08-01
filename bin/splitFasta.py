"""
Splits a fasta file into one file per sequence, and reformats into
a single line. Name of the output file is the name of the sequence.
All files are created in a specified output directory.
E.g. splits this:
    >id1 blah blah
    aactgagtacatgtagacta
    actagtacatagtacgat
    >id2 blah blah
    gagtacatgtagactattgc
    tacatagtac
into two files:
   id1:
    aactgagtacatgtagactaactagtacatagtacgat

   id2:
    gagtacatgtagactattgctacatagtac

This script reads from stdin and writes to files in the specified output directory.
"""
import sys
import os
import argparse
import re

def getOpts () :
  parser = argparse.ArgumentParser()
  parser.add_argument("-c",
    metavar="REGEX",
    dest="chrRegex",
    default=".*",
    help="Regular expression the chromosomes must match.")
  parser.add_argument(
    "-o",
    dest="odir",
    default=".",
    help="Output directory.")
  opts = parser.parse_args()
  opts.chrRegex = re.compile("^(%s)$" % opts.chrRegex)
  return opts

def split (ifd, odir):
  line = ifd.readline()
  ofd = None
  writing = False
  while line:
      if line.startswith('>'):
          seqid = line.split()[0][1:]
          writing = opts.chrRegex.match(seqid)
          if writing:
              ofile = os.path.join(odir, seqid)
              if ofd: ofd.close()
              ofd = open(ofile, 'w')
          line = ifd.readline()
      else:
          if writing:
              ofd.write(line[:-1])
          line = ifd.readline()
  if ofd: ofd.close()

def main ():
  global opts
  opts = getOpts()
  split(sys.stdin, opts.odir)

main()
