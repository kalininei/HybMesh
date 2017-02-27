import command
import hybmeshpack.basic.proc as bp
import hybmeshpack.basic.interf as interf
from hybmeshpack.gdata.framework import Framework


class CommandFlow(bp.AbstractSender):
    # messages from command flow to _send_message procedure
    APPEND_COMMAND = 0
    BEFORE_EXECUTION = 1
    SUCCESS_EXECUTION = 2
    FAILED_EXECUTION = 3
    UNDO_COMMAND = 4
    TO_ZERO_STATE = 5

    def __init__(self):
        super(CommandFlow, self).__init__()
        # command list
        self._commands = [command.Start()]
        # commands[i > curposition] can be redone
        self._curpos = 0
        # commands[i <= startposition] cannot be undone
        self._startpos = 0
        # receiver: gdata.Framework.framework
        self.receiver = Framework()
        # default communication with user
        self.interface = interf.BasicInterface()

    def set_interface(self, interface):
        self.interface = interface
        self.interface.set_flow(self)

    def set_receiver(self, rec):
        self.receiver = rec

    def com_count(self):
        "->int. Number of commands"
        return len(self._commands)

    def com(self, i):
        "-> Command. Get a command by index"
        return self._commands[i]

    # can undo/redo from current position
    def can_undo(self):
        return self._curpos > self._startpos

    def can_redo(self):
        return len(self._commands) - 1 > self._curpos

    # removes all commands after curpos, adds the command to list
    # and executes it
    def exec_command(self, c):
        # remove all commands from current position to last command
        if (self.can_redo()):
            for c in self._commands[self._curpos + 1:]:
                c.reset()
            self._commands = self._commands[:self._curpos + 1]
        # addition
        self.append_command(c)
        # execution
        self.exec_next()

    # Adds the command to the end of the commands list.
    # Doesn't execute it
    def append_command(self, c):
        # append a command
        self._commands.append(c)
        # link current flow to the command
        c.set_flow(self)
        self._send_message(self.APPEND_COMMAND)

    # commands execution procedures
    def exec_next(self):
        "executes (or redo's) the next command"
        if (self.can_redo()):
            # execution
            self._curpos += 1
            res = False

            cmd = self._commands[self._curpos]
            self._send_message(self.BEFORE_EXECUTION, com=cmd)
            try:
                res = cmd.do(self.receiver)
            except Exception as e:
                self.interface.error_handler().execution_error(e, cmd)

            if res:
                self._send_message(self.SUCCESS_EXECUTION, com=cmd)
            else:
                self._curpos -= 1
                self._send_message(self.FAILED_EXECUTION, com=cmd)

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
            self._send_message(self.UNDO_COMMAND)

    def undo_all(self):
        while (self.can_undo()):
            self.undo_prev()

    def purge_flow(self):
        """ resets all commands, deletes commands after current one,
            start position = current position.
        """
        for c in self._commands:
            c.reset()
        self._commands = self._commands[:self._curpos + 1]
        self._startpos = self._curpos

    # sets command receiver to its ZeroState. Commands are preserved
    def to_zero_state(self):
        self._curpos = 0
        self._startpos = 0
        for c in self._commands:
            c.reset()
        self.receiver.to_zero_state()
        self._send_message(self.TO_ZERO_STATE)

    def make_checkpoint(self):
        "returns a checkpoint copy of the current command flow"
        ret = CommandFlow(self._collection)
        ret.receiver = self.receiver.deepcopy()
        ret._curpos = self._curpos
        ret._startpos = self._curpos
        factory = self._collection._factory
        for c in self._commands[1:]:
            ret.append_command(factory.create_from_string(str(c)))
        for i, c in enumerate(self._commands):
            ret.com(i).set_comment(c.get_comment())
        return ret


class FlowCollection(object):
    'Represents the collection of all working flows'
    ACTUAL_FLOW_CHANGED = 0  # Flow newact - new actual flow
    WARNING = 1  # str txt -- warning text
    ERROR = 2  # str txt -- warning text

    def __init__(self):
        """
            Initialises flow collection.
        """
        # functions of arguments (int tp, FlowCollection sender, **kwargs)
        self._subscribers = []
        self._flows = bp.NamedList([("Flow1", CommandFlow(self))])
        _, self._act_flow = self._flows.get_by_index(0)

    def get_actual_flow(self):
        ' -> CommandFlow '
        return self._act_flow

    def set_actual_flow(self, flow_name):
        ' sets flow with flow_name as active '
        self._actFlow = self._flows[flow_name]
        self._actFlow.receiver.view_update()
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
                self._send_message(self.WARNING,
                                   text="Can not remove the last flow")
                return

        # if deleted flow is actual one change actual flow
        if self._flows[flow_name] == self._actFlow:
            ind = self._flows.keys().index(flow_name)
            if ind == 0:
                ind = 1
            else:
                ind -= 1
            self.set_actual_flow(self._flows.keys()[ind])

        # removing from list
        del self._flows[flow_name]

    def rename_flow(self, oldname, newname):
        self._flows.change_key(oldname, newname)

    # function which are called from CommandFlow
    def _new_receiver(self):
        ' returns new empty unconnected framework '
        return self._rec_cls()

    def checkpoint_current(self, flow_name):
        """
            Makes a copy of active flow.
            flowName is the name of new flow
        """
        new_flow = self._actFlow.make_checkpoint()
        self._flows[flow_name] = new_flow
