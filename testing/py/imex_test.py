import os.path
from hybmeshpack import hmscript as hm
hm.check_compatibility("0.2.1")


def check(cond):
    import traceback
    if not cond:
        print "TEST FAILED <<<<<<<<<<<<<<<<<<<<<<<"
        traceback.print_stack()


print "import from fluent msh file"
g1 = hm.import_grid_msh("input1.msh")
check(hm.info_grid(g1) ==
      {'cell_types': {3: 182, 4: 1004},
       'Nnodes': 1196, 'Nedges': 2382, 'Ncells': 1186})

print "export to fluent msh file"
hm.export_grid_msh(g1, "imexout1.msh")
check(os.path.getsize('imexout1.msh') == os.path.getsize('output1.msh'))

print "export to gmsh file"
hm.export_grid_gmsh(g1, "imexout2.msh")
check(os.path.getsize('imexout2.msh') == os.path.getsize('output2.msh'))


print "saving project"
hm.save_project("imexout3.hmp")
