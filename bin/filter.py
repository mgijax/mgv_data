import sys
from lib.Filter import getFilter

for line in getFilter(sys.argv[1], sys.stdin):
    sys.stdout.write(line)
