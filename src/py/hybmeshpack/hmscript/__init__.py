from hybmeshpack.com.flow import CommandFlow


# global flow and data
flow = CommandFlow()


class ExecError(Exception):
    """ Exception that is raised if hmscript function failed """
    def __init__(self, _funname):
        """ _funname - hmscript function name"""
        mes = ''.join(["Error at ", _funname])
        super(ExecError, self).__init__(mes)


class ExportError(Exception):
    """ Exception that is raised if exporting function failed """
    def __init__(self, _ainfo=None):
        """ _ainfo (str) - additional information"""
        mes = _ainfo
        super(ExportError, self).__init__(mes)

class InvalidArgument(ValueError):
    """Raised when user passes invalid argument to a function"""
    def __init__(self, msg=""):
        super(InvalidArgument, self).__init__(msg)

from generalfun import *  # NOQA

from exportfun import *  # NOQA
from importfun import *  # NOQA

from o2info import *  # NOQA
from c2construct import *   # NOQA
from c2oper import *  # NOQA
from g2construct import *  # NOQA
from g2oper import *  # NOQA

from o3info import *  # NOQA
from o3construct import *  # NOQA
from o3oper import *  # NOQA

import _dbg  # NOQA
