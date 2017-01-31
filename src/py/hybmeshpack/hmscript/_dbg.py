import traceback
import sys
from hybmeshpack import progdata
import ctypes as ct
from hybmeshpack.hmcore.proc import cport, ccall


def check(cond):
    if not cond:
        print "TEST FAILED <<<<<<<<<<<<<<<<<<<<<<<"
        traceback.print_stack()


def checkdict(data, kv):
    for k, v in kv.iteritems():
        if k not in data or v != data[k]:
            check(False)


def check_ascii_file(test_hash, fname, opt=""):
    """ Checks ascii file for the hash value.
        all float numbers are rounded, '\r' sybmols deleted before
        hash calculating.
        Possible opt values:
          "" - checks always
          "dev" - checks only if development version is executed
          "linux" - only if unix system
    """
    if opt == "dev":
        if str(progdata.HybMeshVersion.current) != '9.9.9':
            return
    if opt == "linux":
        if not sys.platform.startswith('linux'):
            return
    fname = fname.encode('utf-8')
    ret = ct.c_size_t
    ccall(cport.get_ascii_file_hash, fname, ret)
    ret = ret.value
    check(ret == test_hash)
    if ret != test_hash:
        print "\tinput file hash:", test_hash
        print "\treal file hash: ", ret
