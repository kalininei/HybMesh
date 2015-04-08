""" Describes command flow management """
import copy
import xml.etree.ElementTree as ET
import bp
import globvars


class CommandReceiver(object):
    ' Interface for Command Receiver object '

    def to_zero_state(self):
        'set receiver to its initial state '
        raise NotImplementedError

    def post_proc(self):
        'procedure is called after each Command.do/redo/undo '
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


class Command(object):
    ' abstract base class for flow command '

    def __init__(self, argsdict):
        """ argsdic is a option dictionary {optionname1: optionvalue1, ...}
            each argsdict value should have __str__ method """
        self.__executed = False
        self.__comment = ""
        #comLine is string representation of command:
        #command Code from Factory dictionary + list of arguments
        st = copy.deepcopy(argsdict)
        for k in st.keys():
            st[k] = str(st[k])
        self.__comLine = self._method_code() + ' ' + str(st)
        #last receiver of the command
        self.receiver = None

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
        self.receiver.post_proc()
        return self.__executed

    def undo(self):
        """ undo the command and call receiver.post_proc() method.

            Receiver is taken from last do method call
        """
        if (self.__executed):
            self._undo()
        else:
            raise Exception(str(self) + " can not be reverted")
        self.receiver.post_proc()

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
        return self.__comLine

    #------------ methods to override
    #constructors
    @classmethod
    def fromstring(cls, slist):
        """ -> Command object

            @slist -- string representation of
            constuctor arguments dictionary
        """
        raise NotImplementedError

    #info
    @classmethod
    def _method_code(cls):
        """-> str

           Returns a word which is used to dub this command
           in xml files
        """
        return cls.__name__

    @classmethod
    def doc(cls):
        """-> str.
            Command documentation line as it will be shown in
            commands history window"""
        return cls.__doc__

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


class StartCommand(Command):

    """Start command is the first command in each CommandFlow.

       It reflects the initial state of receiver with no data.
       It should not be executed, redone or undone.
    """

    def __init__(self):
        ' __init__'
        super(StartCommand, self).__init__({})

    def doc(self):
        return "Start"

    @classmethod
    def _method_code(cls):
        return "Start"

    @classmethod
    def fromstring(cls, slist):
        return cls()

    def _exec(self, receiver):
        return True

    def _clear(self):
        pass

    def _undo(self):
        pass

    def _redo(self):
        pass


class Factory(object):
    ' Command class Factory '
    def __init__(self, coms):
        """ Initializes the commands building factory.
            coms - is the list of classes which
            implement Command interface
        """
        #command dictionary: code -> Command class
        self._coms = {}
        #register start command
        self.__register(StartCommand)
        #register user command
        for c in coms:
            self.__register(c)

    def string_from_cls(self, cls):
        ind = self._coms.values().index(cls)
        return self._coms.keys()[ind]

    def cls_from_string(self, s):
        return self._coms[s]

    def create_from_string(self, strcom):
        p = strcom.split(None, 1)
        cls = self.cls_from_string(p[0])
        return cls.fromstring(p[1])

    def create_from_args(self, clscode, *args):
        cls = self.cls_from_string(clscode)
        return cls(*args)

    def __register(self, comcls):
        self._coms[comcls._method_code()] = comcls


class CommandFlow(object):
    #messages from command flow to subscribers
    ExecCommand = 0
    AppendCommand = 1
    UndoCommand = 2
    XmlLoadCommands = 3
    ClearCommands = 4

    def __init__(self, flow_collection):
        #subscribers for command do/undo
        self._subscribers = []
        #reference to collection
        self._collection = flow_collection
        #command list
        self._commands = [StartCommand()]
        #commands[i > curposition] can be redone
        self._curpos = 0
        #commands[i <= startposition] cannot be undone
        self._startpos = 0
        #receiver
        self._receiver = self._collection._new_receiver()

    def com_count(self):
        "->int. Number of commands"
        return len(self._commands)

    def com(self, i):
        "-> Command. Get a command by index"
        return self._commands[i]

    def add_subscriber(self, obj):
        self._subscribers.append(obj)

    def remove_subscriber(self, obj):
        self._subscribers.remove(obj)

    def _subscriber_message(self, tp):
        "sent a message for do/undo subscribers"
        for s in self._subscribers:
            s.command_flow_action(tp)

    #can undo/redo from current position
    def can_undo(self):
        return self._curpos > self._startpos

    def can_redo(self):
        return len(self._commands) > self._curpos + 1

    #removes all commands after curpos, adds the command to list
    #and executes it
    def exec_command(self, c):
        #remove all commands from current position to last command
        if (self.can_redo()):
            self._commands = self._commands[:self._curpos + 1]
        #addition
        self.append_command(c)
        #execution
        self.exec_next()

    #Adds the command to the end of the commands list.
    #Doesn't execute it
    def append_command(self, c):
        self._commands.append(c)
        self._subscriber_message(self.AppendCommand)

    #commands execution procedures
    def exec_next(self):
        "executes (or redo's) the next command"
        if (self.can_redo()):
            #debug mode: save before execution
            if globvars.prog_options.debug_save_before:
                globvars.Flows.xml_save(globvars.prog_options.debug_save_fn)
            #execution
            self._curpos += 1
            if self._commands[self._curpos].do(self._receiver):
                self._subscriber_message(self.ExecCommand)
            else:
                self._curpos -= 1

    def exec_all(self):
        while (self.can_redo()):
            a = self._curpos
            self.exec_next()
            # if no progress was made (error or cancel)
            # stop execution
            if (a == self._curpos):
                break

    def undo_prev(self):
        if (self.can_undo()):
            self._curpos -= 1
            self._commands[self._curpos + 1].undo()
            self._subscriber_message(self.UndoCommand)

    def undo_all(self):
        while (self.can_undo()):
            self.undo_prev()

    #sets command receiver to its ZeroState. Commands are preserved
    def to_zero_state(self):
        self._receiver.to_zero_state()
        self._curpos = 0
        self._startpos = 0
        for c in self._commands:
            c.reset()
        self._receiver.post_proc()
        self._subscriber_message(self.ClearCommands)

    def make_checkpoint(self):
        "returns a checkpoint copy of the current command flow"
        ret = CommandFlow(self._collection)
        ret._receiver = self._receiver.deep_copy()
        ret._curpos = self._curpos
        ret._startpos = self._curpos
        factory = self._collection._factory
        for c in self._commands[1:]:
            ret.append_command(factory.create_from_string(str(c)))
        for i, c in enumerate(self._commands):
            ret.com(i).set_comment(c.get_comment())
        return ret

    def save_state(self, xmlnode):
        #Commands list
        cl = ET.SubElement(xmlnode, "COMMANDS")
        for c in self._commands:
            nd = ET.SubElement(cl, "ENTRY")
            if c is self._commands[self._curpos]:
                nd.attrib['current'] = '1'
            ET.SubElement(nd, "LINE").text = str(c)
            if c.get_comment() != "":
                ET.SubElement(nd, "COMMENT").text = c.get_comment()

        #Program state
        pn = ET.SubElement(xmlnode, "STATE")
        self._receiver.save_state(pn)

    def load_state(self, xmlnode):
        self._receiver.to_zero_state()
        del self._commands[:]
        #Commands list
        cfnd = xmlnode.findall("COMMANDS/ENTRY")
        for cline in cfnd:
            com_text = cline.find("LINE").text
            try:
                com_comment = cline.find("COMMENT").text
            except:
                com_comment = ""
            c = self._collection._factory.create_from_string(com_text)
            c.set_comment(com_comment)
            self.append_command(c)
        self._curpos = cfnd.index(xmlnode.find('COMMANDS/ENTRY[@current]'))
        self._startpos = self._curpos
        #Program state
        pn = xmlnode.find("STATE")
        self._receiver.load_state(pn)
        self._subscriber_message(self.XmlLoadCommands)


class FlowCollection(object):
    ' Represents the collection of all working flows '

    def __init__(self, com_cls, reciever_cls):
        """
            Initialises flow collection.
            @com_cls is the array of classes for each possible command
            @Reciever_cls is the class which implements
            CommandReciever interface
        """
        self._actflow_subscriber = []
        self._factory = Factory(com_cls)
        #set the framework and visualiser class
        self._rec_cls = reciever_cls
        #set flow list and impose active flow
        self._flows = bp.NamedList([("Flow1", CommandFlow(self))])
        _, self._actFlow = self._flows.get_by_index(0)

    def get_actual_flow(self):
        ' -> CommandFlow '
        return self._actFlow

    def set_actual_flow(self, flow_name):
        ' sets flow with flow_name as active '
        self._actFlow = self._flows[flow_name]
        self._actFlow._receiver.view_update()
        for s in self._actflow_subscriber:
            s.actual_flow_changed(flow_name)

    def get_flow(self, name=None, ind=None):
        "-> CommandFlow. Get flow using its name or index"
        if name is not None:
            return self._flows[name]
        elif ind is not None:
            _, ret = self._flows.get_by_index(ind)
            return ret

    def add_actflow_subscriber(self, obj):
        self._actflow_subscriber.append(obj)

    def get_flow_names(self):
        ' -> [All existing flow names]'
        return self._flows.keys()

    def get_actual_flow_name(self):
        ' -> Actual Flow name '
        _, ret = self._flows.get_by_value(self.get_actual_flow())
        return ret

    def remove_flow(self, flow_name):
        ' removes the flow with specified name '
        if len(self._flows) < 2:
                raise Exception("Can not remove the last flow")

        #if deleted flow is actual one change actual flow
        if self._flows[flow_name] == self._actFlow:
            ind = self._flows.keys().index(flow_name)
            if ind == 0:
                ind = 1
            else:
                ind -= 1
            self.set_actual_flow(self._flows.keys()[ind])

        #removing from list
        del self._flows[flow_name]

    def rename_flow(self, oldname, newname):
        self._flows.change_key(oldname, newname)

    def xml_load(self, fname):
        ' fills object from xml file '
        xmlnode = ET.parse(fname).getroot()
        self.__load_states(xmlnode)

    def xml_save(self, fname):
        ' saves all flows to a xml file '
        outp = ET.Element("ComGridProject")
        self.__save_states(outp)
        bp.xmlindent(outp)
        tree = ET.ElementTree(outp)
        tree.write(fname, xml_declaration=True, encoding='utf-8')

    #function which are called from CommandFlow
    def _new_receiver(self):
        ' returns new empty unconnected framework '
        return self._rec_cls()

    def checkpoint_current(self, flow_name):
        """
            Makes a copy of active flow.
            @flowName is the name of new flow
        """
        new_flow = self._actFlow.make_checkpoint()
        self._flows[flow_name] = new_flow

    #private functions
    def __save_states(self, xmlnode):
        ' saves all flows to a xml ElementTree.Element '
        flcol = ET.SubElement(xmlnode, "FLOWS")
        for (nm, f) in self._flows.items():
            fl = ET.SubElement(flcol, "FLOW")
            fl.attrib['name'] = nm
            if f is self._actFlow:
                fl.set("active", "1")
            f.save_state(fl)

    def __load_states(self, xmlnode):
        """ Loads flows from ElementTree.Element.
            All existing flows will be lost.
        """
        #clear all current flows
        self._flows.clear()
        #add new flows
        lst = xmlnode.findall("FLOWS/FLOW")
        #read flow information
        actual_name = lst[0].attrib["name"]
        for fxml in lst:
            nf = CommandFlow(self)
            nf.load_state(fxml)
            nm = fxml.attrib["name"]
            self._flows[nm] = nf
            if ("active" in fxml.attrib):
                actual_name = nm

        self.set_actual_flow(actual_name)
