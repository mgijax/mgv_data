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

def getArgs (cmdLineTokens=None) :
  if not cmdLineTokens: cmdlineTokens = sys.argv
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
  opts = parser.parse_args(cmdLineTokens)
  opts.chrRegex = re.compile("^(%s)$" % opts.chrRegex)
  return opts

def log(s):
    sys.stderr.write(s)

def split (ifd, odir, chr_re, logFcn):
  ofd = None
  writing = False
  lcount = 0
  try:
      line = next(ifd)
      while line:
          lcount += 1
          if line.startswith('>'):
              seqid = line.split()[0][1:]
              writing = chr_re.match(seqid)
              if writing:
                  # output is plain text, not Fasta
                  ofile = os.path.join(odir, seqid + '.txt')
                  if ofd: ofd.close()
                  ofd = open(ofile, 'w')
                  logFcn(line)
              else:
                  logFcn("Skipping: " + line)
              line = next(ifd)
          else:
              if writing:
                  ofd.write(line[:-1])
              line = next(ifd)
  except StopIteration:
      pass
  if ofd: ofd.close()

def main ():
  global opts
  opts = getArgs()
  split(sys.stdin, opts.odir, opts.chrRegex, log)

if __name__ == "__main__":
    main()
