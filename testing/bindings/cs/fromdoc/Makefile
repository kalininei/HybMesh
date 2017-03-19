INCLUDES = -r:System.Windows.Forms.dll -r:System.Drawing.dll

all:
	mcs -debug app.cs ../../../../build/bindings/cs/Hybmesh.cs $(INCLUDES)
	mono --debug app.exe
