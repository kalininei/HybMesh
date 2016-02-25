from hybmeshpack import com, gdata
import hybmeshpack.com.flow  # NOQA


# global flow and data
flow = com.flow.CommandFlow()
data = gdata.Framework()
flow.set_receiver(data)


class ExecError(Exception):
    """ Exception that is raised if hmscript function failed """
    def __init__(self, _funname):
        """ _funname - hmscript function name"""
        mes = ''.join(["Error at ", _funname])
        super(ExecError, self).__init__(mes)


from proto import *   # NOQA
from inout import *  # NOQA
from gproc import *  # NOQA
from cproc import *  # NOQA
from oper import *  # NOQA
from info import *  # NOQA
