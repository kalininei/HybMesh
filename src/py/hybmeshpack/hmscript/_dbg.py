import traceback
import sys
from hybmeshpack import progdata
import ctypes as ct
from hybmeshpack.hmcore import libhmcport


def check(cond):
    if not cond:
        print "TEST FAILED <<<<<<<<<<<<<<<<<<<<<<<"
        traceback.print_stack()


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
    c_fname = fname.encode('utf-8')
    libhmcport.get_ascii_file_hash.restype = ct.c_size_t
    r = libhmcport.get_ascii_file_hash(c_fname)
    check(r == test_hash)
    if r != test_hash:
        print "\tinput file hash:", test_hash
        print "\treal file hash: ", r
