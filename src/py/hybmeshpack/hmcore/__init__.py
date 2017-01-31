' procedures for communication with c libraries'
import ctypes as ct
import os
from hybmeshpack import progdata

_libfn = progdata.get_lib_fn('cport')
cport = ct.cdll.LoadLibrary(os.path.abspath(_libfn))


class CppLibError(Exception):
    def __init__(self, mes):
        super(CppLibError, self).__init__("algorithm error: " + mes)
