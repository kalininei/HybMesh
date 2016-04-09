import traceback


def check(cond):
    if not cond:
        print "TEST FAILED <<<<<<<<<<<<<<<<<<<<<<<"
        traceback.print_stack()


def check_ascii_file(test_hash, fname):
    import ctypes as ct
    from hybmeshpack.hmcore import libhmcport
    c_fname = fname.encode('utf-8')
    libhmcport.get_ascii_file_hash.restype = ct.c_size_t
    r = libhmcport.get_ascii_file_hash(c_fname)
    check(r == test_hash)
    if r != test_hash:
        print "\tinput file hash:", test_hash
        print "\treal file hash: ", r
