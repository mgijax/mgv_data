import os
import sys

DDIR=os.environ['DDIR']
HFILE=os.path.join(DDIR, 'homologies.txt')
fd = open(HFILE, 'r')
cid2hid = {}
for line in fd:
  hid, iid = line.strip().split()
  cid2hid[iid]=hid
  
def feature(f):
  attrs = f[8]
  cid = attrs.get('cID', None)
  hid = cid2hid.get(cid, None)
  if hid:
    attrs['hID'] = hid
  return f
