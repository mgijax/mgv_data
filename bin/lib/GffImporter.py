#
#
# importGff3.py
#
# Imports a GFF3 genome annotation file (a la Ensembl) into form needed by MGV.
#
# Chunking:
#    For efficiency, output can be chunked at specified granularities.
#       chunkSize   meaning
#       0           No chunking. Everything in one file. The file's name is the track name.
#                   It is located in the root directory for the genome.
#       1           Chunk by chromosome. Track name = subdirectory in genome root dir.
#                   Then one subdir per chromosome below this.One file per chromosome.
#       > 1         Chunk size. As above, but output per chromosome is chunked by this size block factor.
# 
import types
import os
import sys
import argparse
import json
import time

try:
    from .gff3lite import Gff3Parser, formatLine, parseLine, parsePragmas
except:
    from gff3lite import Gff3Parser, formatLine, parseLine, parsePragmas

TIMESTAMP = time.asctime(time.localtime(time.time()))

#
class Gff3Importer:
  def __init__(self, gffSource, opts):
    self.opts = opts
    self.parser = Gff3Parser(gffSource, returnHeader=True, returnGroups=True, convertDots='.')
    self.datastream = self.parser.iterate()
    #
    self.seenChrs = {}
    self.seenChrOrder = []
    #
    self.outputDir = None
    self.outputFiles = {}
    self.currFileName = None
    self.currFile = None
    # top-level file name and file handle
    self.tlFileName = None
    self.tlFile = None
    #
    self.processHeader()

  def getOpt(self, n):
    v = getattr(self.opts, n, None)
    if v is None:
      return None
    #
    if v.startswith('##'):
      v = self.pragmas[v[2:]]
    return v

  def processHeader (self) :

    self.pragmas = parsePragmas(next(self.datastream))
    self.outputRootDir = os.path.join(self.opts.outputDir, self.opts.genomePath)
    self.outputDir = os.path.join(self.outputRootDir, "models")
    
  def log(self, s):
    sys.stderr.write(s)
    sys.stderr.flush()

  def ensureDirectory (self, d):
    if not os.path.exists(d):
      os.makedirs(d)

  def getFileName(self, track, chr, blk):
    if chr:
      return os.path.join(self.outputDir, track, chr, str(blk) + '.gff3')
    else:
      return os.path.join(self.outputDir, track, str(blk), '.gff3')

  def getFileHandle(self, fname):
    if fname == self.currFileName:
      # currently open file
      return self.currFile
    elif fname in self.outputFiles:
      # different from current but seen before, reopen in append mode
      self.currFile.close()
      self.currFileName = fname
      self.currFile = open(fname, 'a')
      return self.currFile
    else:
      # first time seeing this fname
      self.ensureDirectory(os.path.dirname(fname))
      self.currFile and self.currFile.close()
      self.currFileName = fname
      self.currFile = open(fname, 'w')
      self.outputFiles[fname] = 0

      fn = os.path.basename(fname)
      dn = os.path.basename(os.path.dirname(fname))
      self.log('%s/%s ' % (dn, fn))
      return self.currFile

  def formatFeature(self, f):
    return formatLine(f)

  def writeGrpToBlk(self, grp, track, chr, blk):
    fname = self.getFileName(track, chr, blk)
    fd = self.getFileHandle(fname)
    for f in grp:
      if 'exons' in f[8]:
        fd.write(self.formatFeature(f))
        self.outputFiles[fname] += 1

  def writeGrp (self, track, grp) :
    if len(grp) == 0: return
    chr = grp[0][0]
    if self.opts.chunkSize == 0 :
      self.writeGrpToBlk(grp, track, None, 0)
    elif self.opts.chunkSize == 1:
      self.writeGrpToBlk(grp, track, chr, 0)
    else:
      gStart = min([f[3] for f in grp])
      gEnd = max([f[4] for f in grp])
      startBlk = gStart // self.opts.chunkSize
      endBlk = gEnd // self.opts.chunkSize
      for blk in range(startBlk, endBlk+1):
        self.writeGrpToBlk(grp, track, chr, blk)

  def writeGene(self, f):
    self.tlFile.write(self.formatFeature(f))

  # Processes a gene model group. For top level features, adds a tCount attribute that contains the
  # number of transcripts for that gene. For transcripts, adds an exons attribute whose value is a list 
  # of exon coordinates. Each exon is encoded as 'offset_length', where offset is equal to
  # exon.start - transcript.start, and length is equal to exon.end - exon.start + 1.
  # Also keeps track of all chromosomes seen in the file and the max coordinate seen on each.
  # Writes genes to genes track and transcripts to transcripts track.
  def processGrp(self, grp):
    # Build index from ID to the features that have a Parent attribute with that ID.
    # Keep track of the chromosome in the file, their order, and the max coordinate of features on each.
    pid2kids = {}
    for f in grp:
      seqid = f[0]
      c = self.seenChrs.get(seqid, None)
      if not c:
        c = { 'name': seqid, 'length': 0 }
        self.seenChrs[seqid] = c
        self.seenChrOrder.append(c)
      c['length'] = max(c['length'], f[4])
      for pid in f[8].get('Parent',[]):
        pid2kids.setdefault(pid, []).append(f)
    # For genes, count transcripts and store in tCount attribute.
    # For transcripts, encode exon coordinates and store in exons attribute
    transcripts = []
    for f in grp:
      fid = f[8].get('ID', None)
      if 'Parent' not in f[8]:
        # top level feature
        f[8]['tCount'] = len(pid2kids.get(fid,[]))
        self.writeGene(f)
      elif fid in pid2kids:
        # mid level feature
        # For efficiency of transfer, encode the exons for a transcript into a col9 attribute.
        # For each exon, just need its offset from the start of the transcript and its length.
        # This requires a lot less space (~time) to transfor than representing each exon as a separate
        # feature with full coordinates. NOTE that we lose exon IDs (and potentially other attributes) 
        # by doing this. Maybe worth it, maybe not?
        exons = list([x for x in pid2kids[fid] if x[2] == 'exon'])
        exons.sort(key=lambda e: e[3])
        eexons = f[8].setdefault('exons',[])
        for e in exons:
          offset = e[3] - f[3]
          length = e[4] - e[3] + 1
          eexons.append('%d_%d' % (offset,length))
        #
        cdss = list([x for x in pid2kids[fid] if x[2] == 'CDS'])
        if cdss:
          cdss.sort(key=lambda x: x[3])
          attrs = cdss[0][8]
          cid = attrs['ID']
          prid = attrs.get('protein_id', attrs.get('curie', ''))
          cstart = cdss[0][3]
          cend = cdss[-1][4]
          f[8]['cds'] = '%s|%s|%d|%d' % (cid, prid, cstart, cend)
        #
        transcripts.append(f)
      else:
        # leaf level feature
        pass
    self.writeGrp('transcripts', transcripts)

  def main (self) :
    self.ensureDirectory(self.outputDir)
    self.ensureDirectory(os.path.join(self.outputDir, 'genes'))
    self.tlFileName = os.path.join(self.outputDir, 'genes', '0.gff3')
    self.tlFile = open(self.tlFileName, 'w')
    for i,grp in enumerate(self.datastream):
      toplevel = self.processGrp(grp)
#
def sanitizeName (n):
    return n.replace('/','').lower()

#
# Returns options object from command line. Has these fields:
#       genomePath - name of genome in paths at Ensembl
#       genome - name/label of genome (as shown to user)
#       taxonid - NCBI taxon id
#       timestamp - timestamp
#       outputDir - when the output gets written
#       chunkSize - for chunking transcripts and exons
def getArgs (cmdLineTokens=None) :
  if not cmdLineTokens: cmdlineTokens = sys.argv
  parser = argparse.ArgumentParser(description='Import genome annotations from a GFF3 file.')
  parser.add_argument('-p','--genomePath',
    metavar='NAME',
    dest='genomePath',
    help='Pathname to use with the genome. By default, this is a sanitized version of the genomeName, but you can override that.')
  parser.add_argument('-g','--genomeName',
    metavar='NAME',
    default='##genome-name',
    dest='genome',
    help='Name of the genome. Required. Either specify the name directly (eg, "C57BL/6J") or provide the name of a pragma in the GFF3 header (eg "##genome-name")')
  parser.add_argument('-x','--taxonid',
    metavar='TAXONID',
    default='##taxonid',
    dest='taxonid',
    help='NCBI taxon identifier.')
  parser.add_argument('-T','--timestamp',
    metavar='TIME',
    dest='timestamp',
    default=TIMESTAMP,
    help='Value to use as the timestamp. Default is the current time. Pretty much any string is OK. Specify a timestamp value or a ##pragma name.')
  parser.add_argument('-d','--outputDirectory',
    metavar='DIR',
    dest='outputDir',
    default='./output',
    help='Output directory.')
  parser.add_argument('-k','--chunkSize',
    type=int,
    metavar='CHUNKSIZE',
    dest='chunkSize',
    default=1,
    help='Transcript file chunk size, in bases. 0 = everything in one file (no chunking); 1 = one chunk per chromosome; >1 = chunk by chromosome and by start position / chunkSize.')
  #
  args = parser.parse_args(cmdLineTokens)
  if not args.genomePath:
    args.genomePath = sanitizeName(args.genome)
  return args
#
if __name__ == '__main__':
  args = getArgs()
  Gff3Importer(sys.stdin, args).main()

