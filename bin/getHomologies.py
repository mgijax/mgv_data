#
# processHomologies.py
#
# Reads JSON files containing homologous ID pairs.
# Puts them all into a single graph of gene IDs and edges.
# Then enumerates the connected components.
# Outputs each gene ID with its connected component ID
#
# usage:
#     python processHomologies.py file1.json file2.json ...
#
import sys
import json

class CCfinder:
  def __init__(self):
    self.neighbors = {}
    self.ccs = []
    if len(sys.argv) == 1:
      self.sources = [sys.argv]
    else:
      self.sources = map(open, sys.argv[1:])

  # Reads the ID pairs in an input file.
  def pairs(self, src):
    j = json.loads(src.read())
    for r in j['results']:
      g1 = r['gene']['id']
      s1 = r['gene']['symbol']
      o1 = r['gene']['taxonId']
      #
      g2 = r['homologGene']['id']
      s2 = r['homologGene']['symbol']
      o2 = r['homologGene']['taxonId']
      #
      yield (g1,g2)

  # Builds the graphs
  def buildGraph (self) :
    for s in self.sources:
      for (a,b) in self.pairs(s):
        self.neighbors.setdefault(a,[]).append(b)
        self.neighbors.setdefault(b,[]).append(a)

  # Recursive traversal step.
  # Args:
  #   n - the current node (a gene ID)
  #   ccc - current connected component (list of nodes)
  def reach(self, n, ccc):
    if n in self.visited:
      return
    self.visited.add(n)
    ccc.add(n)
    for nn in self.neighbors.get(n, []):
      self.reach(nn, ccc)
  #
  def enumerateCCs (self):
    self.visited = set()
    self.count = 0
    self.ccs = []
    for n in self.neighbors.keys():
      ccc = set()
      self.reach(n, ccc)
      if len(ccc): self.ccs.append(ccc)

  #
  def outputCCs (self) :
    for (i,cc) in enumerate(self.ccs):
      for x in cc:
        print (i, x)
  #
  def main(self):
    self.buildGraph()
    self.enumerateCCs()
    self.outputCCs()

#
CCfinder().main()
