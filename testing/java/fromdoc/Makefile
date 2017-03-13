JAVAC = javac
JAVA = java
HMJAVAPATH = /home/ek/HybMesh/build/bindings/java/
HMLIBPATH = /home/ek/HybMesh/build/bin

all:
	$(JAVAC) App.java GridData.java $(HMJAVAPATH)/Hybmesh.java
	$(JAVA) -Djava.library.path=$(HMLIBPATH) -cp ".:$(HMJAVAPATH)" App

