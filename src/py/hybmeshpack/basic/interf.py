import cb  # NOQA


class DefaultErrorHandler(object):
    """ basic implementation of program warnings and errors
        All ErrorHandlers should be inherited from this
    """

    def __init__(self):
        pass

    def execution_error(self, exc, cmd):
        """ called by command flow if command fails to execute
            exc - Exception
            cmd - command which failed
        """
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

    def ask_for_callback(self):
        return cb.SilentCallbackCancel2()
