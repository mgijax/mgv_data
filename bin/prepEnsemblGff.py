#
# prepEnsemblGff.py
#
# Performs conversions on the GFF3 as downloaded from Ensembl in preparation for import into MGV data store.
#

import types
import re
import os
import sys 
from gff3lite import parseLine, formatLine
import argparse
import json
import time

def getOpts () :
  parser = argparse.ArgumentParser()
  parser.add_argument("-c",
    metavar="REGEX",
    dest="chrRegex",
    default=".*",
    help="Regular expression the chromosomes must match.")
  parser.add_argument("-x",
    metavar="TYPE,TYPE,...",
    dest="exclude",
    default="",
    help="SO types to exclude. Matches against values in column 3.")
  opts = parser.parse_args()
  if opts.exclude:
    opts.exclude = opts.exclude.split(',')
  else:
    opts.exclude = []
  opts.chrRegex = re.compile("^(%s)$" % opts.chrRegex)
  return opts

lastLine = None
def write (line) :
  global lastLine
  sys.stdout.write(line)
  lastLine = line

def handleComment (line):
  if line.strip() == "###" and line == lastLine:
    return
  if line.startswith('##sequence-region'):
    chrom = line.split()[1]
    if not opts.chrRegex.match(chrom):
      return
  write(line)

def handleFeature (f) :
  if f[2] in opts.exclude:
    return
  if not opts.chrRegex.match(f[0]):
    return
  '''
  attrs = f[8]
  if 'ID' in attrs:
    attrs['ID'] = attrs['ID'].split(':')[-1]
  if 'Parent' in attrs:
    attrs['Parent'] = list(map(lambda pid: pid.split(':')[-1], attrs['Parent']))
  '''
  write(formatLine(f))

def main () :
  global opts
  opts = getOpts()
  for line in sys.stdin:
    if line.startswith('#'):
      handleComment(line)
    else:
      f = parseLine(line)
      handleFeature (f)

main()
