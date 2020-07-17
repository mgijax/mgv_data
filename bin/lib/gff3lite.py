#
# gff3lite.py
#
# A small GFF3 parser
#
# Reference: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#
# Usage:
#    import Gff3Parser from gff3lite
#    # iterate over the features. Each feature is returned as a list of 9 items.
#    for r in Gff3Parser("myfile.gff3").iterate():
#        # zero-based (GFF spec is 1-based)
#        chr = r[0]
#        # column 9 is a regular dict
#        id = r[8]['ID'] 
#        # coordinates converted to int
#        length = r[4] - r[3] + 1
#
import sys
import types
import urllib.parse as ulib
import json
import re

GFF3HEADER = '##gff-version 3\n'
TAB = '\t'
NL = '\n'
SEMI = ';'
COMMA = ','
EQ = '='
HASH = '#'
GT = '>'

MULTIVALUED = ["Parent", "Dbxref", "Alias", "Note", "Ontology_term"]
GROUPSEPARATOR = "###"

class Gff3Parser :
  def __init__(self, source, returnHeader=True, returnGroups=False, convertDots='.'):
    self.source = source
    self.returnHeader = returnHeader
    self.returnGroups = returnGroups
    self.convertDots = convertDots

  def open(self):
    if type(self.source) is str:
      self.sfd = open(source, 'r')
    else:
      self.sfd = self.source
    self.currLine = None

  def close(self):
    if type(self.source) is str:
      sfd.close()

  def iterate (self):
    self.open()
    header = []
    group = []
    currSeqid = None
    inSeqs = False
    for line in self.sfd:
      if inSeqs:
          continue
      if line.startswith(GT):
          inSeqs = True
          continue
      # Comment line
      if line.startswith(HASH):
        if self.returnGroups and line.strip() == GROUPSEPARATOR:
          # User wants groups and line is a "###" separator
          if group: yield group
          group = []
        elif header != None:
          header.append(line)
        continue
      # Feature line
      if header != None and self.returnHeader:
        # User wants header and this is the first feature. 
        # Yield the header then remember that we did so.
        yield ''.join(header)
        header = None
      #
      record = list([self.convertDots if a == '.' else a for a in parseLine(line)])
      if self.returnGroups:
        if record[0] != currSeqid or "Parent" not in record[8]:
          # if chromosome changes or the feature is top-level (has no parent), start a new group
          if group: yield group
          group = []
        group.append(record)
      else:
        yield record
      currSeqid = record[0]
    # end for loop
    if self.returnGroups and group:
      yield group
    #
    self.close()

  def sortIterate (self) :
    # 
    origRH = self.returnHeader
    origRG = self.returnGroups
    self.returnHeader = True
    self.returnGroups = True
    gffStream = self.iterate()
    header = next(gffStream)
    if origRH:
        yield header
    allModels = list(gffStream)
    # sort by models by chromosome, then start position of the top-level feature
    allModels.sort(key = lambda m: (m[0][0], m[0][3]))
    for model in allModels:
        if origRG:
            yield model
        else:
            for m in model:
                yield m
#
def parseColumn9 (text) :
  c9 = {}
  if text == ".":
    return c9
  parts = text.split(SEMI)
  for p in parts:
    p = p.strip()
    if len(p) == 0:
      continue
    bits = p.split(EQ, 1)
    if len(bits) == 1:
        # syntax error - unescaped ';' in an attribute value?
        # ignore it, and try to keep going
        continue
    n = bits[0].strip()
    v = bits[1].strip()
    if n in MULTIVALUED:
      c9[n] = list(map(ulib.unquote, v.split(COMMA)))
    else:
      c9[n] = ulib.unquote(v)
  return c9

#
def formatColumn9(c9):
  parts = []
  for (n,v) in list(c9.items()):
    nn = ulib.quote(n)
    if type(v) is str:
      vv = ulib.quote(v)
    elif type(v) is list:
      vv = COMMA.join(map(ulib.quote, v))
    else:
      vv = ulib.quote(str(v))
    parts.append("%s=%s" % (nn, vv))
  return SEMI.join(parts)
    
#
def parseLine (line) :
  try:
      if line.endswith(NL) :
        line = line[:-1]
      flds = line.split(TAB)
      if len(flds) != 9:
        raise RuntimeError("Line does not have 9 columns.")
      flds[3] = int(flds[3])
      flds[4] = int(flds[4])
      flds[8] = parseColumn9(flds[8])
      return flds
  except Exception as e:
      sys.stderr.write("Error parsing GFF line: " + line)
      raise

#
PRAGMA_RE = re.compile(r'#[#!]([-\w]+) (.*)')
def parsePragmas(lines):
  ps = {}
  for l in lines:
    m = PRAGMA_RE.match(l)
    if m:
      n = m.group(1)
      v = m.group(2)
      if n in ps:
        if type(ps[n]) is str:
          ps[n] = [ps[n], v]
        else:
          ps[n].append(v)
      else:
        ps[n] = v
  return ps
#
#
def formatLine (row):
  r = row[:]
  r[8] = formatColumn9(r[8])
  return TAB.join(map(str, r)) + NL

if __name__ == "__main__":
  for r in Gff3Parser(sys.stdin,returnGroups=True, returnHeader=True).iterate():
    print(r)
    print("###")

