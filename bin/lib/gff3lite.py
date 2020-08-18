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

  # In nature, one gene can be inside another, e.g., gene A might exist in an intron of gene B.
  # Some providers of GFF3 files segregate each model into its own contiguous set of lines, while
  # others will interleave the lines for gene A within the lines of gene B. Downstream code
  # assumes segregated, so this routine exists to convert an interleaved stream of features
  # into a (segregated) stream of models (lists of features).
  #
  # If the data are already segregated, using this method is safe.
  #
  # Each call passes the next feature in the input (or None to terminate, see below), and returns 
  # a list of zero or more models that have been flushed from the internal buffer and should be output.
  # After the last feature in the input stream, this method should be called one more time
  # with f == None to return the any models remaining in the buffer.
  #
  def deInterleaveNext(self, f):
    #
    if f is None:
        # signals the end of the input. Return anything left in the buffer/
        ret = self.buffer
        self.buffer = []
        self.id2group = {}
        return ret
    #
    attrs = f[8]
    if 'Parent' not in attrs:
        # top level feature.
        # 1. flush items from buffer. Must be careful to preserve ordering of top level features
        flushed = []
        for grp in self.buffer:
            tlf = grp[0]
            # Flush a group if the current feature shows we've moved beyond its top level feature or
            # moved to a different chromosome.
            if tlf[0] != f[0] or tlf[4] < f[3]:
                for g in grp:
                    self.id2group.pop(g[8].get('ID',None), None)
                flushed.append(grp)
            else:
                # don't check the rest of the buffer, or we lose
                # input ordering
                break
        # remove flushed items from the buffer
        self.buffer = self.buffer[len(flushed):]
        # start new group and add to buffer
        grp = [f]
        self.buffer.append(grp)
        if 'ID' in attrs:
            self.id2group[attrs['ID']] = grp
        return flushed
    else:
        # subfeature 
        for p in attrs['Parent']:
            grp = self.id2group.get(p, None)
            if not grp:
                self.pending.setdefault(p,[]).append(f)
            else:
                grp.append(f)
                if 'ID' in attrs:
                    self.id2group[attrs['ID']] = grp
                    grp += self.pending.pop(attrs['ID'],[])
        return []

  #
  def iterate (self) :
    self.open()
    header = []
    group = []
    currSeqid = None
    self.buffer = []
    self.pending = {}
    self.id2group = {}
    for line in self.sfd:
      # detect start of sequence section. Exit loop if found.
      if line.startswith(GT):
          break
      # Comment line
      if line.startswith(HASH):
        if header != None and line.strip() != GROUPSEPARATOR:
          header.append(line)
        continue
      # Non-header line. If this is the first time, yield the header
      if header != None:
        # User wants header and this is the first feature. 
        # Yield the header then remember that we did so.
        yield ''.join(header)
        header = None
      # Feature line. Parse it.
      f = list([self.convertDots if a == '.' else a for a in parseLine(line)])
      # Add feature to buffer. If anything gets flushed, yield it.
      flushed = self.deInterleaveNext(f)
      for grp in flushed:
        yield grp
    # End of input. Yield whatever is left in the buffer/
    for grp in self.deInterleaveNext(None):
        yield grp
    #
    if len(self.pending) > 0 :
        sys.stderr.write("Orphan records detected. " + str(self.pending) + "\n")

  #
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
    if type(r) is str:
        print(r)
    else:
        for f in r:
            print(formatLine(f), end = '')
        print("###")

