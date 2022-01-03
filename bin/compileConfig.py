#
# compileConfig.py
#
import sys
import os
import yaml
import json
from lib.UrlResolver import resolveAll

with open(sys.argv[1], "r") as file:
    config = yaml.safe_load(file)

MY_DIR=os.path.dirname(os.path.realpath(__file__))
BCORE = os.path.join(MY_DIR, "buildcore.sh")

print("#!/usr/bin/bash")
print("source %s" % BCORE)
print("activateVenv")
print("parseCommandLine $*")
print('logit "Command line: $CMDLINE"')
for gc in config["buildList"]:
    resolveAll(gc)
    print("#")
    print("# %s" % gc["name"])
    print("#")
    print('if [[ %s == $MATCHPAT ]] ; then' % gc["path"])
    print("\texport GCONFIG='%s'" % json.dumps(gc))
    print("\tGDIR='%s'" % gc["path"])
    for dc in gc["tracks"]:
        print('\tif [[ $TRACK == "%s" || $TRACK == "all" ]] ; then' % dc["type"])
        print("\t\tTTYPE='%s'" % dc["type"])
        print("\t\tFURL='%s'" % dc["url"])
        print("\t\tFTYPE='%s'" % dc["filetype"])
        print("\t\texport DCONFIG='%s'" % json.dumps(dc))
        if gc["type"] == "homology" :
            fcmd = 'python %s/importAllianceOrthology.py' % MY_DIR
            print("\t\tFILTER='%s'" % fcmd)
            print("\t\timportHomology")
        elif gc["type"] == "genome":
            fcmd = 'cat' if not 'filter' in dc else ('python %s/filter.py %s' % (MY_DIR, dc["filter"]))
            print("\t\tFILTER='%s'" % fcmd)
            print("\t\timportGenome")
        else:
            raise RuntimeError("Unknown type: " + gc["type"])
        print("\tfi")
    print("fi")
