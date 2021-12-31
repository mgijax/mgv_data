import sys
import importlib
import json
import os

GCONFIG = json.loads(os.environ.get("GCONFIG", "{}"))
DCONFIG = json.loads(os.environ.get("DCONFIG", "{}"))
sys.stderr.write('GCONFIG=' + str(GCONFIG) + '\n')
sys.stderr.write('DCONFIG=' + str(DCONFIG) + '\n')

farg = sys.argv[1]
clsname = farg[0].capitalize() + farg[1:]
modname = "lib.%s" % clsname
mod = importlib.import_module(modname)
fcls = getattr(mod, clsname)
filt = fcls(sys.stdin, GCONFIG, DCONFIG)

for line in filt:
    sys.stdout.write(line)
