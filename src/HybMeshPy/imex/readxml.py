'module imlements all read from xml routines'
# import xml.etree.ElementTree as ET


# ------------ Commands
def load_command_flow(comflow, xmlnode):
    'adds commands to flow.CommandFlow comflow from xmlnode tagged FLOW'
    from com import factory
    #Commands list
    cfnd = xmlnode.findall("COMMANDS/ENTRY")
    for cline in cfnd:
        com_title = cline.find("TITLE").text
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
    curcom = xmlnode.find('COMMANDS/ENTRY[@current]')
    if curcom is not None:
        comflow._curpos = cfnd.index(curcom)


# ------------- Frameworks
def load_framework_state(fw, xmlnode):
    'adds data to framework from xmlnode tagged STATE'
    import gdata
    #boundary types
    bnode = xmlnode.find("BTYPES")
    if bnode is not None:
        fw.boundary_types.xml_load(bnode)

    #Grids
    gnodes = xmlnode.findall("GRIDS/GRID2")
    for nd in gnodes:
        name = nd.attrib["name"]
        grd = gdata.grid2.Grid2.create_from_xml(nd)
        fw.add_grid(name, grd)

    #User contours
    cnodes = xmlnode.findall("CONTOURS/CONTOUR2")
    for nd in cnodes:
        name = nd.attrib['name']
        cont = gdata.contour2.AbstractContour2.create_from_xml(nd)
        fw.add_user_contour(name, cont)
