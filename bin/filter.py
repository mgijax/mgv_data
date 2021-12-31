import sys
import importlib

farg = sys.argv[1]
clsname = farg[0].capitalize() + farg[1:] + "Filter"
modname = "lib.%s" % clsname
mod = importlib.import_module(modname)
fcls = getattr(mod, clsname)
filt = fcls(sys.stdin)

for line in filt:
    sys.stdout.write(line)
