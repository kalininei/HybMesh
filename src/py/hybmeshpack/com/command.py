""" Describes abstract command interface """
import ast
import copy
from hybmeshpack import basic
import hybmeshpack.basic.proc as bp


class CommandReceiver(object):
    ' Interface for Command Receiver object '

    def to_zero_state(self):
        'set receiver to its initial state '
        raise NotImplementedError

    def save_state(self, xmlnode):
        'save to xml node'
        raise NotImplementedError

    def load_state(self, xmlnode):
        'load from xml node'
        raise NotImplementedError

    def deep_copy(self):
        'make a deep copy of current receiver state'
        raise NotImplementedError


class ExecutionError(bp.EmbException):
    def __init__(self, message, sender=None, upper_error=None):
        '(text message, Command object)'
        super(ExecutionError, self).__init__(message, upper_error)
        self.message = message
        self.sender = sender

    def __str__(self):
        s = 'HybMesh error in command\n'
        if self.sender is not None:
            s += str(self.sender) + '\n'
        else:
            s += 'Unknown\n'
        s += '\nMessage:\n'
        s += self.message
        return s


class ObjectNotFound(ExecutionError):
    def __init__(self, objname, sender=None, upper_error=None):
        s = "Object %s not found" % objname
        super(ObjectNotFound, self).__init__(s, sender, upper_error)


class BasicOption(object):
    """ Represent option transformation before
        str(Command.option) will be called.
        For types like int, float, str, [int], {int: int}, ...
        Base class will be enough.
        For complicated ones child class is needed.
        See example of child implementation for gdata.Point2 option below.
    """
    def __init__(self, tp=None):
        self.tp = tp

    def serial(self, v):
        'before writing'
        return v

    def unserial(self, s):
        'after reading'
        if self.tp is None:
            return s
        else:
            return self.tp(s)


class ListOfOptions(BasicOption):
    "list of simple data: int, floats, str"
    def __init__(self, tp):
        'tp - another BasicOption object'
        super(ListOfOptions, self).__init__(tp)

    def serial(self, v):
        return [self.tp.serial(x) for x in v]

    def unserial(self, v):
        return [self.tp.unserial(x) for x in v]


class NoneOr(BasicOption):
    "simple data or None"
    def __init__(self, tp):
        'tp - another BasicOption object'
        super(NoneOr, self).__init__(tp)

    def serial(self, v):
        return None if v is None else self.tp.serial(v)

    def unserial(self, v):
        if v is None or v == 'None':
            return None
        else:
            return self.tp.unserial(v)


class SubDictOption(BasicOption):
    def __init__(self, **kwargs):
        """ Present dictionary of options:
            name1: value1,
            name2: value2,...

            kwargs are  name=BasicOption()
        """
        self.args = kwargs

    def serial(self, v):
        a = {}
        for k, val in v.items():
            if k in self.args:
                a[k] = self.args[k].serial(val)
        return a

    def unserial(self, v):
        a = {}
        for k, val in v.items():
            if k in self.args:
                a[k] = self.args[k].unserial(val)
        return a


class Point2Option(BasicOption):
    def __init__(self):
        from hybmeshpack.basic.geom import Point2
        super(Point2Option, self).__init__(Point2)

    def serial(self, v):
        return (v.x, v.y)

    def unserial(self, v):
        return self.tp(v[0], v[1])


class BoolOption(BasicOption):
    def __init__(self):
        super(BoolOption, self).__init__(bool)

    def serial(self, v):
        return v

    def unserial(self, v):
        if v in [0, False, '0', 'False', 'false']:
            return False
        else:
            return True


class ListCompressedOption(BasicOption):
    def __init__(self):
        super(ListCompressedOption, self).__init__()

    def serial(self, v):
        return bp.compress_int_list(v)

    def unserial(self, v):
        return bp.int_list_from_compress(v)


class Command(object):
    ' abstract base class for flow command '

    def __init__(self, argsdict):
        """ argsdic is a option dictionary {optionname1: optionvalue1, ...}
            Its to-string, from-string methods are defined in _arguments_types
            child class overloads
        """
        self.__executed = False
        self.__comment = ""
        self.options = copy.deepcopy(argsdict)

        #last receiver of the command
        self.receiver = None

        # flow to which this command belongs
        # It is used for getting interface callback only
        self.parent_flow = None

    def set_flow(self, flow):
        self.parent_flow = flow

    def do(self, receiver):
        """ execute (or redos) the command
            returns True if success
        """
        if (self.__executed):
            self._redo()
        else:
            self.receiver = receiver
            if (self._exec()):
                self.__executed = True
        return self.__executed

    def undo(self):
        """ undo the command

            Receiver is taken from last do method call
        """
        if (self.__executed):
            self._undo()
        else:
            raise Exception(str(self) + " can not be reverted")

    def reset(self):
        ' clears all backups. Undo operation is not possible after reset call '
        if (self.__executed):
            self.__executed = False
            self._clear()

    #user comments property
    def set_comment(self, s):
        ' sets string comment to the command '
        self.__comment = s

    def get_comment(self):
        ' returns string comment of the command '
        return self.__comment

    def __str__(self):
        return self.method_code() + " | " + self.opt_line()

    def opt_line(self):
        '->str. Used for xml writing of the command options'
        a = {}
        rules = self._arguments_types()
        for k, v in self.options.items():
            a[k] = rules[k].serial(v)
        return str(a)[1:-1]

    @classmethod
    def method_code(cls):
        """-> str

           Returns a word which is used to dub this command
           in xml files
        """
        return cls.__name__

    #constructors
    @classmethod
    def fromstring(cls, slist):
        """ -> Command object

            slist -- string representation of
            constuctor arguments dictionary
        """
        a = ast.literal_eval('{' + slist + '}')
        tps = cls._arguments_types()
        for (k, v) in tps.items():
            if k in a and a[k] is not None:
                a[k] = v.unserial(a[k])
        return cls(a)

    #------------ methods to override
    #info
    @classmethod
    def doc(cls):
        """-> str.
            Command documentation line as it will be shown in
            commands history window
        """
        return cls.__doc__

    @classmethod
    def _arguments_types(cls):
        """ returns dictionary of {argname: CommandOption}
            for simple types like int, float, ... which could be
            serialized whitout string converting to-string and from-string
            data are not
        """
        raise NotImplementedError

    #evalution
    def _exec(self):
        """ command execution returns True
            if success and False otherwise
        """
        raise NotImplementedError

    def _clear(self):
        raise NotImplementedError

    def _undo(self):
        raise NotImplementedError

    def _redo(self):
        raise NotImplementedError

    #auxilliary functions
    def cont_by_name(self, name):
        '->contour by its name or raise ObjNotFound'
        try:
            _, _, c = self.receiver.get_ucontour(name=name)
            return c
        except:
            raise ObjectNotFound(name, self)

    def grid_by_name(self, name):
        '->grid by its name or raise ObjNotFound'
        try:
            _, _, c = self.receiver.get_grid(name=name)
            return c
        except:
            raise ObjectNotFound(name, self)

    def any_by_name(self, name):
        '->(grid, contour) or raise. One of (grid, contour) is None'
        if name in self.receiver.get_grid_names():
            return self.receiver.get_grid(name=name)
        elif name in self.receiver.get_ucontour_names():
            return self.receiver.get_ucontour(name=name)
        else:
            raise ObjectNotFound(name, self)

    def any_cont_by_name(self, name):
        '->contour or raise. Grid contour or user contour by object name'
        if name in self.receiver.get_grid_names():
            return self.receiver.get_grid(name=name)[2].cont
        elif name in self.receiver.get_ucontour_names():
            return self.receiver.get_ucontour(name=name)[2]
        else:
            raise ObjectNotFound(name, self)

    def ask_for_callback(self, tp=None):
        """ ask parent flow interface for callback.
            Returns proper callback object
            tp - callback types from interf.
               = None -> CB_CANCEL2 callback
        """
        try:
            return self.parent_flow.get_interface().ask_for_callback(tp)
        except:
            return basic.interf.Callback.silent_factory(tp)


class Start(Command):
    """Start command is the first command in each CommandFlow.
       It reflects the initial state of receiver with no data.
       It should not be executed, redone or undone.
    """

    def __init__(self):
        ' __init__'
        super(Start, self).__init__({})

    def doc(self):
        return "Start"

    @classmethod
    def _arguments_types(cls):
        return {}

    def _exec(self, receiver):
        return True

    def _clear(self):
        pass

    def _undo(self):
        pass

    def _redo(self):
        pass
