# Filters out biological_regions and chromosomes
# Gets mouse/rat orthologs (as MGIid/RGDid pairs) from MouseMine.
# For each gene, extracts RGD id (if it exists), uses it to find MGI id of mouse ortholog (if it exists)
# and sets canonical id to that MGI id.
# Converts column 3 to 'protein_coding_gene' where appropriate.

import sys
import gff3lite
import urllib
import re
import argparse

#
HUMANMINE="http://www.humanmine.org/humanmine/service"
MOUSEMINE="http://www.mousemine.org/mousemine/service"
#
HGNC_re = re.compile(r'Source:HGNC.*Acc:(HGNC:\d+)')
RGD_re = re.compile(r'Source:RGD.*Acc:(\d+)')
#
def getOpts () :
  parser = argparse.ArgumentParser(description="Prepare human or rat genome file from Ensembl.")
  parser.add_argument('-g', '--organism',
    dest="organism",
    required=True,
    choices=["human","rat"],
    help="The organism. One of human, rat.")
  return parser.parse_args()
#
def doMineQuery(wsurl, q):
    fmt = 'tab'
    url = '%s/query/results?format=%s&query=%s' % (wsurl,fmt,urllib.quote_plus(q))
    fd = urllib.urlopen(url)
    for line in fd:
        toks = line[:-1].split('\t')
        yield toks
    fd.close()
#
# Given a query that returns two columns (ie, pairs of IDs), runs the query and returns an index
# from the first columns' value to the second.
def doIndexQuery(wsurl, q, unique=True, reversed=False, mapper=None):
  index = {}
  lhsIndex = 0 if reversed else 1
  rhsIndex = 1 - lhsIndex
  for r in doMineQuery(wsurl, q):
    if mapper:
      r = mapper(r)
    if unique:
      index[r[lhsIndex]] = r[rhsIndex]
    else:
      index.setdefault(r[lhsIndex],[]).append(r[rhsIndex])
  return index
#
# Returns an index from Ensembl human gene ID to MGI id of mouse ortholog
# Have to get this from HumanMine
def getMgiEnsembl () :
  q = '''<query model="genomic" view="Gene.primaryIdentifier Gene.homologues.homologue.secondaryIdentifier">
    <constraint path="Gene.organism.taxonId" code="A" op="=" value="10090"/>
    <constraint path="Gene.homologues.homologue.organism.taxonId" code="B" op="=" value="9606"/>
    <constraint path="Gene.primaryIdentifier" code="C" op="CONTAINS" value="MGI:"/>
    <constraint path="Gene.homologues.type" code="D" op="==" value="least diverged orthologue"/>
  </query>
  '''
  return doIndexQuery(HUMANMINE, q)
#
# Returns an index from HGNC human gene ID to MGI id of mouse ortholog
def getMgiHgnc () :
  q = '''
  <query model="genomic"
    view="Gene.primaryIdentifier Gene.homologues.homologue.crossReferences.identifier">
    <constraint path="Gene.organism.taxonId" code="A" op="=" value="10090"/>
    <constraint path="Gene.homologues.homologue.organism.taxonId" code="B" op="=" value="9606"/>
    <constraint path="Gene.homologues.homologue.crossReferences.source.name" code="C" op="=" value="HGNC"/>
    <constraint path="Gene.homologues.dataSets.name" code="E" op="=" value="Mouse/Human Orthologies from MGI"/>
  </query>
  '''
  return doIndexQuery(MOUSEMINE, q)
#
# Returns an index from RGD rat gene ID to MGI id of mouse ortholog
def getMgiRgd () :
  q = '''
  <query model="genomic"
    view="Gene.primaryIdentifier Gene.homologues.homologue.primaryIdentifier">
    <constraint path="Gene.organism.taxonId" code="A" op="=" value="10090"/>
    <constraint path="Gene.homologues.homologue.organism.taxonId" code="B" op="=" value="10116"/>
    <constraint path="Gene.homologues.type" code="C" op="==" value="least diverged orthologue"/>
  </query>
  '''
  return doIndexQuery(MOUSEMINE, q, mapper=lambda r: [r[0], r[1].replace('RGD:','')])
#
def processRecord (rec, index, regexp):
  if rec[2] in ['biological_region', 'chromosome']:
    return None
  attrs = rec[8]
  attrs2 = {}
  #
  ensemblid = attrs.get('ID', None)
  if ensemblid:
    ensemblid = attrs2['ID'] = ensemblid.split(':')[-1] # remove the 'GeneID:' prefix
  #
  parents = attrs.get('Parent', None)
  if parents:
    attrs2['Parent'] = map(lambda p: p.split(':')[-1], parents)
  #
  name = attrs.get('Name', None)
  if name:
    attrs2['symbol'] = name
  #
  if 'biotype' in attrs:
    bt = attrs2['biotype'] = attrs['biotype']
    if bt == 'protein_coding':
      rec[2] = 'protein_coding_gene'
  # description=DEAD/H-box helicase 11 like 1 [Source:HGNC Symbol%3BAcc:HGNC:37102];
  description = attrs.get('description', '')
  m = regexp.search(description)
  xid = m.group(1) if m else None
  if ensemblid:
    cid = index.get(ensemblid, None)
    if not cid:
      cid = index.get(xid, None)
    if cid:
      attrs2['cID'] = cid

  rec[8] = attrs2
  return rec

def main () :
  opts = getOpts()
  if opts.organism == "human":
    index = getMgiEnsembl()
    index.update(getMgiHgnc())
    regexp = HGNC_re
  else:
    index = getMgiRgd()
    regexp = RGD_re
  for line in sys.stdin:
    if line.startswith('#'):
      sys.stdout.write(line)
      continue
    rec = processRecord(gff3lite.parseLine(line), index, regexp)
    if rec:
      sys.stdout.write(gff3lite.formatLine(rec))

#
main()

