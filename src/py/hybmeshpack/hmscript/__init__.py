import traceback
import os
import sys
import struct
from hybmeshpack.com.flow import CommandFlow
import hybmeshpack.basic.interf
import hybmeshpack.basic.cb
from hybmeshpack.basic.cb import UserInterrupt
import hybmeshpack.basic.proc
from decorator import decorator


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


@decorator
def hmscriptfun(method, *args, **kw):
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
            newe = ExecError(method.__name__, str(e))
            raise newe.__class__, newe, exc_info[2]
    return ret(*args, **kw)


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
        super(ConsoleErrorHandler, self).execution_error(exc, cmd)


class ConsoleInterface0(hybmeshpack.basic.interf.BasicInterface):
    """ silent mode interface """
    def __init__(self):
        super(ConsoleInterface0, self).__init__()

    # overriden functions
    def flow_messages(self, tp, sender, **kwargs):
        if tp == sender.FAILED_EXECUTION:
            # exception will be catched in hmscript function
            raise Exception
        elif tp == sender.USER_INTERRUPT:
            # exception will be catched in hmscript function
            raise UserInterrupt


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
        self.show_linenum = True

    # overriden functions
    def flow_messages(self, tp, sender, **kwargs):
        if tp == sender.BEFORE_EXECUTION:
            line = kwargs['com'].method_code()
            if self.show_linenum:
                for s in traceback.extract_stack()[::-1]:
                    if not s[0].startswith('<') and "hybmeshpack" not in s[0]:
                        fn = os.path.basename(s[0])
                        line = "{} ({}:{})".format(line, fn, s[1])
                        break
            print line
        elif tp == sender.SUCCESS_EXECUTION:
            print '--- OK'
        elif tp == sender.USER_INTERRUPT:
            print '--- Interrupted'
        elif tp == sender.FAILED_EXECUTION:
            print '--- Failed'
        super(ConsoleInterface2, self).flow_messages(
            tp, sender, **kwargs)


class ConsoleInterface3(ConsoleInterface2):
    """ verbosity=3 mode interface """
    def __init__(self):
        super(ConsoleInterface3, self).__init__()

    def ask_for_callback(self):
        return hybmeshpack.basic.cb.ConsoleCallbackCancel2()


def console_interface_factory(verbosity):
    if verbosity == 0:
        return ConsoleInterface0
    elif verbosity == 1:
        return ConsoleInterface1
    elif verbosity == 2:
        return ConsoleInterface2
    elif verbosity == 3:
        return ConsoleInterface3
    else:
        raise Exception("Imposible vertosity: {}".format(verbosity))


class PipeInterface(hybmeshpack.basic.interf.BasicInterface):

    def __init__(self, pipe_read, pipe_write, verbosity=0):
        super(PipeInterface, self).__init__()
        self.p = [pipe_read, pipe_write]
        self.last_error_message = ""
        self.messenger = console_interface_factory(verbosity)()
        self.messenger.show_linenum = False
        self.callback_base = self.messenger.ask_for_callback().__class__
        self.error_handler_base = self.messenger.error_handler().__class__

    def reset(self, verbosity):
        return PipeInterface(*self.p, verbosity=verbosity)

    def ask_for_callback(self):
        class LCallback(self.callback_base):
            def __init__(self, parent):
                super(LCallback, self).__init__()
                self.parent = parent

            def _callback(self, n1, n2, p1, p2):
                l1, l2 = len(n1), len(n2)
                ilen = 8 + 8 + 4 + 4 + l1 + l2
                s = struct.pack('=iddii%is%is' % (l1, l2),
                                ilen, p1, p2, l1, l2, n1, n2)
                os.write(self.parent.p[1], "B")
                os.write(self.parent.p[1], s)
                r = os.read(self.parent.p[0], 1)
                if r == "S":
                    self._proceed = False
                elif r == "G":
                    self._proceed = True
                return super(LCallback, self)._callback(n1, n2, p1, p2)

        return LCallback(self)

    # overriden functions
    def error_handler(self):

        class LErrorHandler(self.error_handler_base):
            def __init__(self, parent):
                super(LErrorHandler, self).__init__()
                self.parent = parent

            def execution_error(self, exc, cmd):
                self.parent.last_error_message = str(exc)
                return super(LErrorHandler, self).execution_error(exc, cmd)

        return LErrorHandler(self)

    # overriden functions
    def flow_messages(self, tp, sender, **kwargs):
        if tp == sender.FAILED_EXECUTION:
            # exception will be catched in hmscript function
            raise Exception(self.last_error_message)

        self.messenger.flow_messages(tp, sender, **kwargs)


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
