import ctypes as ct


class UserInterrupt(Exception):
    """Raised when callback returns 1"""
    def __init__(self):
        super(UserInterrupt, self).__init__("Interrupted by user")


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
        if tp is int:
            return ct.c_int
        elif tp is str:
            return ct.c_char_p
        elif tp is float:
            return ct.c_double

    def __construct_c_function(self, py_func):
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

    # this should be called from py side procedures
    # if callback object is availible.
    def pycall(self, *args):
        if self._callback(*args) == 1:
            raise UserInterrupt()

    # custruct callback delegate for py or c functions.
    # Its return should be handled. To raise immediately use pycall
    def _get_callback(self):
        def cb(*args):
            return self._callback(*args)
        if self._is_c_target:
            return self.__construct_c_function(cb)
        else:
            return cb

    # function for overloading
    def _callback(self, *args):
        return 0 if self._proceed else 1

    def subcallback(self, part, total):
        return SubCallback(part, total, self)

    @staticmethod
    def silent_factory():
        return SilentCallbackCancel2()


class SubCallback(Callback):
    def __init__(self, part, total, main):
        super(SubCallback, self).__init__(main.cb_res_type, main.cb_args_types)
        self.__main = main
        self.__part = part
        self.__total = total

    def _callback(self, n1, n2, p1, p2):
        p1 = float(self.__part + p1) / self.__total
        return self.__main._callback(n1, n2, p1, p2)


class SilentCallbackCancel2(Callback):
    """ Default implementation of callback function of type
        Callback.CB_CANCEL2
        Provides no callback output
    """
    def __init__(self):
        super(SilentCallbackCancel2, self).__init__(
            int, (str, str, float, float))


class ConsoleCallbackCancel2(SilentCallbackCancel2):
    """ Callback to console. no cancel.
    """
    def __init__(self):
        super(ConsoleCallbackCancel2, self).__init__()
        self.__prev_n1, self.__prev_n2 = '', ''

    def subcallback(self, part, total):
        self.__prev_n1, self.__prev_n2 = '', ''
        return super(ConsoleCallbackCancel2, self).subcallback(part, total)

    @staticmethod
    def _supl(s, char, sz):
        if len(s) >= sz:
            return s
        add = char * (sz - len(s))
        return s + add

    @staticmethod
    def _breakline(s, sz):
        ret = []
        icur = 0
        while icur < len(s):
            ret.append(s[icur: icur + sz])
            icur += sz
        return ret

    def _callback(self, n1, n2, p1, p2):
        # [--n--] -s- [--n--]
        n, s = 30, 4
        par1, par2 = '', ''
        len1, len2 = int(n * p1), int(n * p2)

        if n1 and n1 != self.__prev_n1:
            par1 = n1
            self.__prev_n1 = n1
        if n2 and n2 != self.__prev_n2:
            par2 = n2
            self.__prev_n2 = n2

        if len(par1) > n:
            par1 = par1[:n-3] + '...'
        if len(par2) > n:
            par2 = par2[:n-3] + '...'
        par1 = self._supl(par1, '#', len1)
        par1 = self._supl(par1, '-', n)
        if n2 != "Done":
            par2 = self._supl(par2, '#', len2)
        par2 = self._supl(par2, '-', n)
        ret = ''.join(['[', par1, ']', s * ' ', '[', par2, ']'])
        print ret
        return super(ConsoleCallbackCancel2, self)._callback(n1, n2, p1, p2)
