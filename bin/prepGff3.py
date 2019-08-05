#
# prepGff3.py
#
# Performs conversions on the GFF3 as downloaded in preparation for import into MGV data store.
# Universal operations like filtering by type are performed here. Other specific changes are made
# via modules specified on the command line (-m). Each module exports a function named feature()
# that takes one feature as argument, makes any desired changes, then returns it, or None. Returning
# None causes that feature to be skipped.
#

import types
import re
import os
import sys 
from gff3lite import parseLine, formatLine
import argparse
import json
import time
import importlib

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
  parser.add_argument("-m",
    metavar="NAME,NAME...",
    dest="modules",
    help='''Names of modules (comma separated, no spaces) to use for munging 
      each GFF3 record.
      Each module defines a callable named "feature" which takes a feature
      as argument, does watever, and returns a feature or None. 
      Returning None causes the feature to be skipped.
      Modules are run in the order listed.''')
  opts = parser.parse_args()
  #
  if opts.exclude:
    opts.exclude = opts.exclude.split(',')
  else:
    opts.exclude = []
  #
  opts.chrRegex = re.compile("^(%s)$" % opts.chrRegex)
  #
  if opts.modules:
    opts.modules = list(map(lambda m: importlib.import_module(m), opts.modules.split(',')))
  else:
    opts.modules = []
  #
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
  for m in opts.modules:
    if not m.feature(f):
      return
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
