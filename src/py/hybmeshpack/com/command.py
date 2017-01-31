""" Describes abstract command interface """
import ast
import copy
from hybmeshpack import basic
import hybmeshpack.basic.proc as bp


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


class Command(object):
    ' abstract base class for flow command '

    def __init__(self, argsdict):
        """ argsdic is a option dictionary {optionname1: optionvalue1, ...}
            Its structure is defined in _arguments_types
            child class overloads
        """
        self.__executed = False
        self.__comment = ""
        self.options = copy.deepcopy(argsdict)

        #last receiver of the command of gdata.Framework type
        self.receiver = None

        # flow to which this command belongs (flow.Commandflow)
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
        if self.__executed:
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
           in xml files and error outputs
        """
        return cls.__name__

    #constructors
    @classmethod
    def fromstring(cls, slist):
        """ -> Command object

            slist -- string representation of
            constuctor arguments dictionary
        """
        if slist is None:
            slist = ""
        a = ast.literal_eval('{' + slist + '}')
        tps = cls._arguments_types()
        for (k, v) in tps.items():
            if k in a:
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
    def get_option(self, code):
        try:
            return self.options[code]
        except KeyError:
            if self._arguments_types[code].has_default():
                return self._arguments_types[code].default
            else:
                raise KeyError("Option %s is uknown" % code)

    def cont2_by_name(self, name):
        '->contour by its name or raise ObjNotFound'
        try:
            return self.receiver.get_contour2(name)
        except KeyError:
            raise ObjectNotFound(name, self)

    def grid2_by_name(self, name):
        '->grid by its name or raise ObjNotFound'
        try:
            return self.receiver.get_grid2(name)
        except KeyError:
            raise ObjectNotFound(name, self)

    def grid3_by_name(self, name):
        '->grid by its name or raise ObjNotFound'
        try:
            return self.receiver.get_grid3(name)
        except KeyError:
            raise ObjectNotFound(name, self)

    def surface3_by_name(self, name):
        '->surface by its name or raise ObjNotFound'
        try:
            return self.receiver.get_surface3(name)
        except KeyError:
            raise ObjectNotFound(name, self)

    def any_cont_by_name(self, name):
        """ ->Contour2 or raise. Grid contour or user contour by object name.
        Makes a deep copy if name is a grid.
        """
        return self.any_acont_by_name().contour2()

    def any_acont_by_name(self, name):
        """ ->AbstractContour2 or raise.
        """
        try:
            return self.receiver.any_contour(name)
        except KeyError:
            raise ObjectNotFound(name, self)

    def any_asurface_by_name(self, name):
        """ ->AbstractSurface3 or raise.
        """
        try:
            return self.receiver.any_surface(name)
        except KeyError:
            raise ObjectNotFound(name, self)

    def any_surface_by_name(self, name):
        """ ->Surface3 or raise. Grid surface or surface by object name
        Makes a deep copy if name is a grid.
        """
        return self.any_asurface_by_name().surface3()

    def ask_for_callback(self):
        """ ask parent flow interface for callback.
            Returns proper callback object
        """
        try:
            return self.parent_flow.get_interface().ask_for_callback()
        except:
            return basic.interf.Callback.silent_factory()


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
