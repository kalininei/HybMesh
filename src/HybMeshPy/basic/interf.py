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


class SilentCallbackCancel2PG(object):
    """ Default implementation of callback function of type

        int callback_fun(str BaseName, str SubName,
                double proc1, double proc2)
        it returns 1 for cancellation requiry and 0 otherwise
        BaseName/SubName - names of main/sub process
        proc1/proc2 - [0,1] progress of base/sub processes

        No callback, No output here
    """
    def __init__(self):
        pass

    def initialize(self, func, args, **kwargs):
        """ func(*args, callback_fun) - target procedure,
        """
        self._func = func
        self._args = args + (self._get_callback(),)
        self._result = None
        self._proceed = True

    def execute_command(self):
        self._result = None
        self._result = self._func(*self._args)

    def get_result(self):
        return self._result

    def _get_callback(self):
        import ctypes as ct

        def cb(n1, n2, p1, p2):
            self._info(n1, n2, p1, p2)
            return 0 if self._proceed else 1

        #if target function is a c function then convert callback to
        #a c function pointer
        if isinstance(self._func, ct._CFuncPtr):
            cbfunc = ct.CFUNCTYPE(ct.c_int, ct.c_char_p, ct.c_char_p,
                   ct.c_double, ct.c_double)
            cb2 = cbfunc(cb)
            return cb2
        else:
            return cb

    def _info(self, n1, n2, p1, p2):
        pass


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

    def ok_cancel_2pg_callback(self):
        return SilentCallbackCancel2PG()


class ConsoleCallbackCancel2PG(SilentCallbackCancel2PG):
    """
        Callback to console.
    """
    def __init__(self):
        super(ConsoleCallbackCancel2PG, self).__init__()
        self.__prev_n1, self.__prev_n2 = '', ''

    def _info(self, n1, n2, p1, p2):
        n, s = 25, 4
        if n1 != self.__prev_n1 or n2 != self.__prev_n2:
            outs = ''
            if n1 != self.__prev_n1:
                outs = n1
            if len(outs) < n + s + 2:
                outs += ' ' * (n + s + 2 - len(outs))
            if n2 != self.__prev_n2:
                outs += n2
            self.__prev_n1 = n1
            self.__prev_n2 = n2
            print outs

        progress1, progress2 = ['-'] * n, ['-'] * n
        w1, w2 = int(n * p1), int(n * p2)
        for i in range(w1):
            progress1[i] = '#'
        for i in range(w2):
            progress2[i] = '#'
        progress1 = ['['] + progress1 + [']']
        progress2 = ['['] + progress2 + [']']

        print ''.join(progress1) + ' ' * s + ''.join(progress2)


class ConsoleInterface(BasicInterface):
    'Verbose console interface'

    def __init__(self):
        super(ConsoleInterface, self).__init__()

    def flow_messages(self, tp, sender, **kwargs):
        'self.cur_flow messages receiver'
        if tp == sender.BEFORE_EXECUTION:
            print kwargs['com'].method_code()
        if tp == sender.SUCCESS_EXECUTION:
            print '--- OK'
        if tp == sender.FAILED_EXECUTION:
            print '--- Failed'

    def flow_collection_messages(self, tp, sender, **kwargs):
        'self.cur_flow_collection messages receiver'
        pass

    def ok_cancel_2pg_callback(self):
        return ConsoleCallbackCancel2PG()
