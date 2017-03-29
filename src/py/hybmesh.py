#!/usr/bin/env python2
"""
Console interface for HybMesh. Possible arguments:
-v -- print version and exit
-u -- check for updates and exit
-sx fn.py [-silent] [-verbosity n]
    Execute hybmesh python script fn.py
    -silent: no console callback during execution
    -verbosity n: choose console output mode
        0: silent mode (same as -silent)
        1: +execution errors descriptions
        2: +commands start/end reports
        3: +progress bar [default]
-px p1 p2 [-silent]
    Execute hybmesh in a pipe communication mode.
    p* are open pipe descriptors for reading (p1) and writing (p2)
"""

# -x fn.hmp [-sgrid gname fmt fn] [-sproj fn] [-silent]
#     Execute command flow from 'hmp' file and saves resulting data.
#     Options are
#     -sgrid gname fmt fn: export resulting grid with called gname to
#     file fn using format fmt. Can be called multiple times.
#         Formats:
#             vtk - vtk format
#             hmg - native HybMesh format
#             msh - fluent mesh format

#     -sproj fn.hmp: save project file after execution to file fn.hmp
#     -silent: no console callback during execution

#     Example:
#     > hybmesh -x state1.hmp -sgrid Grid1 vtk grid1.vtk -sproj final.hmp
#         1) loads data from stat1.hmp
#         2) executes command flow
#         3) exports Grid1 to vtk format
#         4) saves finalized project to final.hmp
import sys
from hybmeshpack import progdata, imex, basic


def sxexec(argv):
    from hybmeshpack import hmscript
    fn = sys.argv[argv.index('-sx') + 1]
    verb = 3
    if '-silent' in argv:
        verb = 0
    elif '-verbosity' in argv:
        try:
            iv = argv.index('-verbosity')
            verb = int(argv[iv + 1])
            if verb < 0 or verb > 3:
                raise
        except:
            sys.exit('Invalid verbosity level. See -help.')
    hmscript.flow.set_interface(hmscript.console_interface_factory(verb)())
    execfile(fn)
    if verb > 1:
        print "DONE"
    sys.exit()


def xexec(argv):
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


def pxexec(argv):
    """ pipe mode execution """
    from hybmeshpack import hmscript
    import os
    import struct
    import ctypes as ct
    pipe_read = int(argv[2])
    pipe_write = int(argv[3])
    if os.name == 'nt':
        import msvcrt
        pipe_read = msvcrt.open_osfhandle(pipe_read, os.O_APPEND | os.O_RDONLY)
        pipe_write = msvcrt.open_osfhandle(pipe_write, os.O_APPEND)
    if '-silent' in argv:
        hmscript.flow.set_interface(hmscript.ConsoleInterface0())
    else:
        hmscript.flow.set_interface(hmscript.PipeInterface(
                pipe_read, pipe_write))
    while 1:
        # wait for a signal from client
        s1 = os.read(pipe_read, 1)
        # command was written to dataread pipe
        if s1 == "C":
            sz = os.read(pipe_read, 4)
            sz = struct.unpack('=i', sz)[0]
            cm = os.read(pipe_read, sz)
            sz = os.read(pipe_read, 4)
            sz = struct.unpack('=i', sz)[0]
            args = os.read(pipe_read, sz)
            try:
                ret = eval("hmscript.{}({})".format(cm, args))
            except hmscript.UserInterrupt:
                # interrupted by callback function
                os.write(pipe_write, "I")
            except Exception as e:
                # error return
                s = str(e)
                os.write(pipe_write, "E")
                os.write(pipe_write, struct.pack('=i', len(s)) + s)
            else:
                os.write(pipe_write, "R")
                # normal return
                if ret is None:
                    # None return
                    os.write(pipe_write, '\0\0\0\0')
                elif isinstance(ret, ct.Array):
                    # print("this is", ret[:])
                    # command returning ctypes array object
                    rawlen, alen = ct.sizeof(ret), ret._length_
                    os.write(pipe_write, struct.pack('=ii', rawlen + 4, alen))
                    os.write(pipe_write, ret)
                else:
                    # regular command which returns python types
                    s = repr(ret)
                    os.write(pipe_write, struct.pack('=i', len(s)) + s)
        # quit signal
        elif s1 == "Q" or not s1:
            sys.exit()
        else:
            raise Exception("Invalid server instruction")


def main():
    if len(sys.argv) == 1 or sys.argv[1] in ['help', 'h', '-help',
                                             '-h', '--help']:
        print __doc__
        sys.exit()

    if len(sys.argv) == 2:
        if sys.argv[1] == '-v':
            print progdata.HybMeshVersion.current()
            sys.exit()
        if sys.argv[1] == '-u':
            r = progdata.check_for_updates()
            print 'Current version: %s' % r[0]
            if r[1] is None:
                print 'Failed to check for latest version'
            elif r[2] == -1:
                print 'New version %s is available at %s' % \
                    (r[1], progdata.project_url())
            else:
                print 'No updates are available'
            sys.exit()
    # execution
    if len(sys.argv) >= 3 and '-sx' in sys.argv:
        sxexec(sys.argv)
    elif len(sys.argv) >= 3 and '-x' in sys.argv:
        xexec(sys.argv)
    elif len(sys.argv) >= 4 and sys.argv[1] == '-px':
        pxexec(sys.argv)

    sys.exit('Invalid options string. See -help.')


if __name__ == '__main__':
    main()
