#!/usr/bin/env python
import subprocess

hm = '../src/HybMeshPy/HybMesh.py'
arg = [
    hm,
    '-x', 'intro_test.hmp',
    '-sproj', '~final_intro_test.hmp',
    '-sgrid', 'GridForExport', 'vtk', '~export1.vtk',
    '-sgrid', 'GridForExport', 'ggen', '~export1.ggen',
    '-sgrid', 'GridForExport', 'msh', '~export1.msh',
    '-sgrid', 'GridForExport', 'hmg', '~export1.hmg',
]

print 'Test Start'
subprocess.call(arg)

