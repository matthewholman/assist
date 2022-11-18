#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob
import subprocess
ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()

with open("version.txt") as f:
    assistversion = f.readlines()[0].strip()
    print("Updating version to "+assistversion)

#with open("doc/index.rst") as f:
#    index = f.readlines()
#
#with open("doc/index.rst","w") as f:
#    for i in range(0,len(index)):
#        if "Welcome to ASSIST" in index[i]:
#            index[i] = "Welcome to ASSIST ("+assistversion+")\n"
#            underline = ""
#            for j in range(len(index[i])-1):
#                underline += "="
#            underline += "\n"
#            index[i+1] = underline
#        f.write(index[i])

with open("README.rst") as f:
    readme = f.readlines()

with open("README.rst","w") as f:
    for i in range(0,len(readme)):
        if "![Version]" in readme[i]:
            readme[i] = "[![Version](https://img.shields.io/badge/assist-v"+assistversion+"-green.svg?style=flat)](https://assist.readthedocs.org)\n"
        f.write(readme[i])

with open("src/assist.c") as f:
    assistlines = f.readlines()
    for i,l in enumerate(assistlines):
        if "**VERSIONLINE**" in l:
            assistlines[i] = "const char* assist_version_str = \""+assistversion+"\";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.\n"

    with open("src/assist.c", "w") as f:
        f.writelines(assistlines)

with open("setup.py") as f:
    setuplines = f.readlines()
    for i,l in enumerate(setuplines):
        if "version='" in l:
            setuplines[i] = "    version='"+assistversion+"',\n"
        if "GITHASHAUTOUPDATE" in l:
            setuplines[i] = "    ghash_arg = \"-DASSISTGITHASH="+ghash+"\" #GITHASHAUTOUPDATE\n"

    with open("setup.py", "w") as f:
        f.writelines(setuplines)

#shortversion = assistversion
#while shortversion[-1] != '.':
#    shortversion = shortversion[:-1]
#    
#shortversion = shortversion[:-1]
#
#with open("doc/conf.py") as f:
#    conflines = f.readlines()
#    for i,l  in enumerate(conflines):
#        if "version =" in l:
#            conflines[i] = "version = '"+shortversion+"'\n"
#        if "release =" in l:
#            conflines[i] = "release = '"+assistversion+"'\n"
#
#    with open("doc/conf.py", "w") as f:
#        f.writelines(conflines)
