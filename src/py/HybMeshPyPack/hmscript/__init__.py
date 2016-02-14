from HybMeshPyPack import com, gdata
import HybMeshPyPack.com.flow


# global flow and data
flow = com.flow.CommandFlow()
data = gdata.Framework()
flow.set_receiver(data)

from proto import *   # NOQA
from inout import *  # NOQA
from gproc import *  # NOQA
from cproc import *  # NOQA
from oper import *  # NOQA


