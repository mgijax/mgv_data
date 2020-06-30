#
# getHomologies.py
#
# Processes Alliance orthology download file.
# - Strips out all but the essential bits
# - Splits into one file per taxon id
#
import sys
import os.path
import argparse
import json

taxon2file = {}

#
def getFile (taxonid) :
  global opts
  fname = os.path.join(opts.odir, taxonid + '.json')
  if taxonid not in taxon2file:
    taxon2file[taxonid] = {
      "fname" : fname,
      "fd" : open(fname, 'w'),
      "count" : 0
    }
  return taxon2file[taxonid]
  
#
def closeAll () :
  for rec in list(taxon2file.values()):
    rec["fd"].write("]")
    rec["fd"].close()
#
def outputRecord (id1, tx1, id2, tx2, yn) :
  rec = getFile(tx1)
  if rec["count"] == 0:
      rec["fd"].write("[")
  else:
      rec["fd"].write(",")
  rec["count"] += 1
  rec["fd"].write(json.dumps([id1, tx1, id2, tx2, yn]))
  rec["fd"].write("\n")

#
inCount = 0
def processLine (line) :
  global inCount
  if line.startswith("#"):
    return
  # skip column labels
  inCount += 1
  if inCount == 1:
    return
  # parse and ouput
  fs = line[:-1].split('\t')
  # Example desired record - just the essentials:
  # ['FB:FBgn0033981', 'NCBITaxon:7227', 'MGI:3646373', 'NCBITaxon:10090', 'YN']
  outputRecord(
    fs[0],
    fs[2].replace('NCBITaxon:',''),
    fs[4],
    fs[6].replace('NCBITaxon:',''),
    fs[11][0] + fs[12][0])

#
def getOpts () :
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-d',
    '--directory',
    dest="odir",
    default=".",
    help="Output directory.")
  return parser.parse_args()

#
def main () :
  global opts
  opts = getOpts()
  for line in sys.stdin:
    processLine(line)
  closeAll()

main()
