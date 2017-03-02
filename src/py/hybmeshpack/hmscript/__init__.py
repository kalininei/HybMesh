import traceback
import os
import sys
import struct
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


class UserInterrupt(Exception):
    """Raised when callback returns 1"""
    def __init__(self):
        super(UserInterrupt, self).__init__("Interrupted by user")


def hmscriptfun(method):
    """ decorator for interface functions
    for catching exceptions.
    """
    def ret(*args, **kw):
        try:
            return method(*args, **kw)
        except ExecError:
            raise
        except InvalidArgument:
            raise
        except UserInterrupt:
            raise
        except Exception as e:
            exc_info = sys.exc_info()
            ############################################3
            import traceback
            print traceback.format_exc(exc_info[2])
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
            # exception will be catched in hmscript function
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


class PipeInterface(ConsoleInterface0):
    def __init__(self, sig_read, sig_write, data_read, data_write):
        super(PipeInterface, self).__init__()
        self.p = [sig_read, sig_write, data_read, data_write]
        self.was_cancelled = False
        self.last_error_message = ""

    class Callback(hybmeshpack.basic.cb.SilentCallbackCancel2):
        def __init__(self, interf):
            super(PipeInterface.Callback, self).__init__()
            self.interf = interf
            self.interf.was_cancelled = False

        def _callback(self, n1, n2, p1, p2):
            l1, l2 = len(n1), len(n2)
            ilen = 8 + 8 + 4 + 4 + l1 + l2
            s = struct.pack('=iddii%is%is' % (l1, l2),
                            ilen, p1, p2, l1, l2, n1, n2)
            os.write(self.interf.p[3], s)
            os.write(self.interf.p[1], "B")
            r = os.read(self.interf.p[0], 1)
            if r == "S":
                self._proceed = False
                self.was_cancelled = True
            elif r == "G":
                self._proceed = True
            return super(PipeInterface.Callback, self)._callback()

    class ErrorHandler(hybmeshpack.basic.interf.DefaultErrorHandler):
        def __init__(self, interf):
            super(PipeInterface.ErrorHandler, self).__init__()
            self.interf = interf

        def execution_error(self, exc, cmd):
            #######################
            ConsoleErrorHandler().execution_error(exc, cmd)
            self.interf.last_error_message = str(exc)

    def ask_for_callback(self):
        return self.Callback(self)

    # overriden functions
    def error_handler(self):
        return self.ErrorHandler(self)

    # overriden functions
    def flow_messages(self, tp, sender, **kwargs):
        if tp == sender.FAILED_EXECUTION:
            # exception will be catched in hmscript function
            if self.was_cancelled:
                raise UserInterrupt()
            else:
                raise Exception(self.last_error_message)


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
import _bindoper  # NOQA
import _dbg  # NOQA
