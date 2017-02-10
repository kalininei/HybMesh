import traceback
import os
import sys
from hybmeshpack.com.flow import CommandFlow
import hybmeshpack.basic.interf
import hybmeshpack.basic.cb
import hybmeshpack.basic.proc


class ExecError(Exception):
    """Raised if hmscript function fails"""
    def __init__(self, funname, msg):
        """ _funname - hmscript function name"""
        mes = "Error at " + repr(funname)
        if msg:
            mes = mes + "\nInfo: " + msg
        super(ExecError, self).__init__(mes)


class InvalidArgument(ValueError):
    """Raised when user passes invalid argument to a hmscript function"""
    def __init__(self, msg=""):
        super(InvalidArgument, self).__init__(msg)


def hmscriptfun(method):
    """ decorator for interface functions
    for catching exceptions.
    """
    def ret(*args, **kw):
        try:
            return method(*args, **kw)
        except InvalidArgument:
            raise
        except Exception as e:
            exc_info = sys.exc_info()
            newe = ExecError(method.__name__, str(e))
            raise newe.__class__, newe, exc_info[2]
    return ret


class ConsoleErrorHandler(hybmeshpack.basic.interf.DefaultErrorHandler):
    def __init__(self):
        super(ConsoleErrorHandler, self).__init__()

    def execution_error(self, exc, cmd):
        print '<<<<<<<<<<<<<<<<<<<<<<<'
        # hybmesh command description
        s = 'HybMesh error in command\n'
        s += str(cmd) + '\n'
        print s
        # traceback
        print ''.join(traceback.format_exc()), '>>>>>>>>>>>>>>>>>>>>>>>'


class ConsoleInterface0(hybmeshpack.basic.interf.BasicInterface):
    """ silent mode interface """
    def __init__(self):
        super(ConsoleInterface0, self).__init__()

    # overriden functions
    def flow_messages(self, tp, sender, **kwargs):
        if tp == sender.FAILED_EXECUTION:
            # exception will be catched in hmscipt function
            raise Exception


class ConsoleInterface1(ConsoleInterface0):
    """ verbosity=1 mode interface """
    def __init__(self):
        super(ConsoleInterface1, self).__init__()

    # overriden functions
    def error_handler(self):
        return ConsoleErrorHandler()


class ConsoleInterface2(ConsoleInterface1):
    """ verbosity=2 mode interface """
    def __init__(self):
        super(ConsoleInterface2, self).__init__()

    # overriden functions
    def flow_messages(self, tp, sender, **kwargs):
        if tp == sender.BEFORE_EXECUTION:
            line = kwargs['com'].method_code()
            for s in traceback.extract_stack()[::-1]:
                if "hybmeshpack" not in s[0]:
                    fn = os.path.basename(s[0])
                    line = "%s (%s:%s)" % (line, fn, s[1])
                    break
            print line
        if tp == sender.SUCCESS_EXECUTION:
            print '--- OK'
        if tp == sender.FAILED_EXECUTION:
            print '--- Failed'
        super(ConsoleInterface2, self).flow_messages(
            tp, sender, **kwargs)


class ConsoleInterface3(ConsoleInterface2):
    """ verbosity=3 mode interface """
    def __init__(self):
        super(ConsoleInterface3, self).__init__()

    def ask_for_callback(self):
        return hybmeshpack.basic.cb.ConsoleCallbackCancel2()


# global flow and data
flow = CommandFlow()

# IMPORTS
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
