#!/usr/bin/env python
"""
Console interface for HybMesh. Possible arguments:
-v: print version and exit
-u: check for updates and exit
-x fn.hmp [-sgrid gname fmt fn] [-sproj fn] [-silent]
    Execute command flow from 'hmp' file and saves resulting data. Options are
    -sgrid gname fmt fn: export resulting grid with called gname to
    file fn using format fmt. Can be called multiple times.
        Formats:
            vtk - vtk format
            hmg - native HybMesh format
            msh - fluent mesh format

    -sproj fn.hmp: save project file after execution to file fn.hmp
    -silent: no console callback during execution

    Example:
    > HybMesh -x state1.hmp -sgrid Grid1 vtk grid1.vtk -sproj final.hmp
        1) loads data from stat1.hmp
        2) executes command flow
        3) exports Grid1 to vtk format
        4) saves finalized project to final.hmp
-sx fn.py [-silent]
    Execute hybmesh python script fn.py
    -silent: no console callback during execution
"""
import sys
from HybMeshPyPack import progdata, imex, basic

def main():
    if len(sys.argv) == 1 or sys.argv[1] in ['help', 'h', '-help', '-h']:
        print __doc__
        sys.exit()

    if len(sys.argv) == 2:
        if sys.argv[1] == '-v':
            print progdata.program_version()
            sys.exit()
        if sys.argv[1] == '-u':
            r = progdata.check_for_updates()
            print 'Current version: %s' % r[0]
            if r[2] is None:
                print 'Failed to check for latest version'
            elif r[2] != 1:
                print 'No updates are availible'
            else:
                print 'Latest version: %s is availible at %s' % (
                            r[1], progdata.project_url())
            sys.exit()

    # execution
    if len(sys.argv) >= 3 and '-sx' in sys.argv:
        fn = sys.argv[sys.argv.index('-sx') + 1]
        from HybMeshPyPack import hmscript
        if '-silent' not in sys.argv:
            hmscript.flow.set_interface(basic.interf.ConsoleInterface())
        execfile(fn)
        if '-silent' not in sys.argv:
            print "DONE"
        sys.exit()
    elif len(sys.argv) >= 3 and '-x' in sys.argv:
        fn = sys.argv[sys.argv.index('-x') + 1]
        # load flow and state
        f = imex.read_flow_and_framework_from_file(fn)

        if '-silent' not in sys.argv:
            f.set_interface(basic.interf.ConsoleInterface())

        # run till end
        f.exec_all()

        # save current state
        if '-sproj' in sys.argv:
            fn = sys.argv[sys.argv.index('-sproj') + 1]
            imex.write_flow_and_framework_to_file(f, fn)

        # export grid
        for i, op in enumerate(sys.argv):
            if op == '-sgrid':
                gname = sys.argv[i + 1]
                fmt = sys.argv[i + 2]
                fn = sys.argv[i + 3]
                try:
                    imex.export_grid(fmt, fn, flow=f, name=gname)
                    print '%s saved to %s' % (gname, fn)
                except Exception as e:
                    print 'Exporting Failure: %s' % str(e)

        print "DONE"
        sys.exit()

    sys.exit('Invalid option string. See -help.')
	


if __name__ == '__main__':
    main()
