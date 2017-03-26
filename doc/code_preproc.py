"""
Script is used to process the source file
before pasting it into documentation.

Usage:
    python code_preproc.py path/to/filename

Creates a resulting file in a current directory called
code_preproc_out-filename

Functionality.

1. Removes lines in block

<code_preproc remove vvv>
....
<code_preproc remove ^^^>
"""
import sys
import os

if len(sys.argv) < 2:
    raise Exception("processed file is not defined")

fn = sys.argv[1]
outfn = "code_preproc_out-" + os.path.split(fn)[1]

fin = open(sys.argv[1], "r")
fout = open(outfn, "w")
canwrite = True

for ln in fin.readlines():
    if "<code_preproc remove vvv>" in ln:
        canwrite = False
    elif "<code_preproc remove ^^^>" in ln:
        canwrite = True
    elif canwrite:
        fout.write(ln)

fin.close()
fout.close()
