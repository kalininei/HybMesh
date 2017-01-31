import traceback
import os
import proc


class DefaultErrorHandler(object):
    """ basic implementation of program warnings and errors
        Flushes everything to console.
        All ErrorHandlers should be inherited from this
    """

    def __init__(self):
        pass

    def known_execution_error(self, ke):
        'ke - com.command.ExecutionError'
        print '<<<<<<<<<<<<<<<<<<<<<<<'
        proc.EmbException.print_estack(ke)
        print ''
        print str(ke)
        print '>>>>>>>>>>>>>>>>>>>>>>>'

    def unknown_execution_error(self, e, com=None):
        'e - Exception'
        print '<<<<<<<<<<<<<<<<<<<<<<<'
        proc.EmbException.print_estack(e)
        if com is not None:
            print ''
            print 'Received from: ' + str(com.method_code())
        print '>>>>>>>>>>>>>>>>>>>>>>>'


class Callback(object):
    """ Abstract callback function of type
        cb_result callback_fun(*cb_args)
    """
    def __init__(self, cb_res_type, cb_args_types):
        self.cb_res_type = cb_res_type
        self.cb_args_types = cb_args_types
        self._proceed = True

    @staticmethod
    def __pytype_to_c_type(tp):
        import ctypes as ct
        if tp is int:
            return ct.c_int
        elif tp is str:
            return ct.c_char_p
        elif tp is float:
            return ct.c_double

    def __construct_c_function(self, py_func):
        import ctypes as ct
        args = ()
        if self.cb_res_type is not None:
            args = args + (self.__pytype_to_c_type(self.cb_res_type),)
        for t in self.cb_args_types:
            args = args + (self.__pytype_to_c_type(t),)
        cbfunc = ct.CFUNCTYPE(*args)
        return cbfunc(py_func)

    def initialize(self, func, args):
        """ func(*args, callback_fun) - target procedure,
        """
        import ctypes as ct
        self._func = func
        self._is_c_target = isinstance(self._func, ct._CFuncPtr)
        self._args = args + (self._get_callback(),)
        self._result = None

    def execute_command(self):
        self._result = None
        self._proceed = True
        self._result = self._func(*self._args)

    def get_result(self):
        return self._result

    def _get_callback(self):
        def cb(*args):
            return self._callback(*args)
        if self._is_c_target:
            return self.__construct_c_function(cb)
        else:
            return cb

    # function for overloading
    def _callback(*args):
        raise NotImplementedError

    def subcallback(self, part, total):
        return self

    @staticmethod
    def silent_factory():
        return SilentCallbackCancel2()


class SilentCallbackCancel2(Callback):
    """ Default implementation of callback function of type
        Callback.CB_CANCEL2
        Provides no callback output
    """
    def __init__(self):
        super(SilentCallbackCancel2, self).__init__(
            int, (str, str, float, float))

    def _callback(self, n1, n2, p1, p2):
        return 0 if self._proceed else 1


class BasicInterface(object):
    'Basic silent interface. Should be parent for all other interfaces'

    def __init__(self):
        self.cur_flow = None
        self.cur_flow_collection = None

    def set_flow(self, flow):
        #unregister from previous flow messages
        if self.cur_flow is not None:
            self.cur_flow.remove_subscriber(self.flow_messages)
        #register for cur_flow messages
        self.cur_flow = flow
        self.cur_flow.add_subscriber(self.flow_messages)

    # ============ Functions for overriding
    def flow_messages(self, tp, sender, **kwargs):
        'self.cur_flow messages receiver'
        pass

    def flow_collection_messages(self):
        'self.cur_flow_collection messages receiver'
        pass

    def error_handler(self):
        return DefaultErrorHandler()

    def ask_for_callback(self):
        return SilentCallbackCancel2()


class ConsoleCallbackCancel2Sub(SilentCallbackCancel2):
    def __init__(self, part, total, main):
        super(ConsoleCallbackCancel2Sub, self).__init__()
        self.__main = main
        self.__part = part
        self.__total = total

    def _callback(self, n1, n2, p1, p2):
        p1 = float(self.__part + p1) / self.__total
        return self.__main._callback(n1, n2, p1, p2)


class ConsoleCallbackCancel2(SilentCallbackCancel2):
    """ Callback to console. no cancel.
    """
    def __init__(self):
        super(ConsoleCallbackCancel2, self).__init__()
        self.__prev_n1, self.__prev_n2 = '', ''

    def subcallback(self, part, total):
        self.__prev_n1, self.__prev_n2 = '', ''
        return ConsoleCallbackCancel2Sub(part, total, self)

    def _callback(self, n1, n2, p1, p2):
        n, s = 25, 4
        if n1 != self.__prev_n1 or n2 != self.__prev_n2:
            outs = ''
            outsnext = []
            if n1 != self.__prev_n1:
                outs = n1
            critlen = n + s + 2
            if len(outs) > critlen:
                fk = outs.split()
                outsnext = [fk[0]]
                for substr in fk[1:]:
                    if len(outsnext[-1] + substr) < critlen:
                        outsnext[-1] += " " + substr
                    else:
                        outsnext.append(substr)
                outs = outsnext[0]
                outsnext = outsnext[1:]
            if len(outs) < critlen:
                outs += ' ' * (critlen - len(outs))
            if n2 != self.__prev_n2:
                outs += n2
            self.__prev_n1 = n1
            self.__prev_n2 = n2
            print outs
            for p in outsnext:
                print p

        progress1, progress2 = ['-'] * n, ['-'] * n
        w1, w2 = int(n * p1), int(n * p2)
        for i in range(w1):
            progress1[i] = '#'
        for i in range(w2):
            progress2[i] = '#'
        progress1 = ['['] + progress1 + [']']
        progress2 = ['['] + progress2 + [']']

        print ''.join(progress1) + ' ' * s + ''.join(progress2)
        return super(ConsoleCallbackCancel2, self)._callback(n1, n2, p1, p2)


class ConsoleInterface(BasicInterface):
    'Verbose console interface'

    def __init__(self):
        super(ConsoleInterface, self).__init__()

    def flow_messages(self, tp, sender, **kwargs):
        'self.cur_flow messages receiver'
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

    def flow_collection_messages(self, tp, sender, **kwargs):
        'self.cur_flow_collection messages receiver'
        pass

    def ask_for_callback(self):
        return ConsoleCallbackCancel2()
