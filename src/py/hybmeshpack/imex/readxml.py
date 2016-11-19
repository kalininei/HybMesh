'module imlements all read from xml routines'
from hybmeshpack.com import factory


# ------------ Commands
def load_command_flow(comflow, xmlnode):
    'adds commands to flow.CommandFlow comflow from xmlnode tagged FLOW'
    #Commands list
    cfnd = xmlnode.findall("COMMANDS/COM")
    # reset all commands because they can not be undone after this procedure
    if len(cfnd) > 0:
        comflow.purge_flow()
    for cline in cfnd:
        try:
            com_title = cline.attrib["name"]
        except:
            continue
        try:
            com_string = cline.find("LINE").text
        except:
            com_string = ""
        try:
            com_comment = cline.find("COMMENT").text
        except:
            com_comment = ""
        if com_title != 'Start':
            c = factory.create_from_string(com_title, com_string)
            comflow.append_command(c)
        else:
            c = comflow._commands[0]
        c.set_comment(com_comment)
        if 'current' in cline.attrib:
            comflow._startpos = len(comflow._commands) - 1
            comflow._curpos = comflow._startpos
        if 'end_load_here' in cline.attrib:
            break


# ------------- Frameworks
def load_framework_state(fw, xmlnode):
    'adds data to framework from xmlnode tagged STATE'
    #boundary types
    if xmlnode is not None:
        fw.boundary_types.xml_load(xmlnode)
