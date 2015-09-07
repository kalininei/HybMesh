#!/usr/bin/env python
import subprocess

hm = '../src/HybMeshPy/HybMesh.py'
arg = [
    hm,
    '-x', 'dbg.hmp',
    '-sproj', '~dbg.hmp',
]

print 'Debug Start'
subprocess.call(arg)
