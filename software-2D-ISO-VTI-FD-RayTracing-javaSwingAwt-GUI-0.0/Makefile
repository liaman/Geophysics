GCC = gcc
JAVAC = javac
JAVAH = javah
JAVA = java

WARNING = -Xlint:deprecation

IJAVAINCLUDE = -I$JAVA_HOME/include
IJAVAINCLUDElinux = -I$JAVA_HOME/include/linux

CFILE = lib/seismic/src_seismic_RungeKutta.c
OLIB = lib/seismic/libRungeKutta.so
LIBPATH = lib/seismic/

MAINCLASS = MainWindow

ALL:
	$(JAVAC) MainWindow.java $(WARNING)
	$(JAVAH) -d ./lib/ src.seismic.RungeKutta
	$(GCC) -shared -fpic -o $(OLIB) $(IJAVAINCLUDE) $(IJAVAINCLUDElinux) $(CFILE)
	$(JAVA) -Djava.library.path=$(LIBPATH) $(MAINCLASS) 
clean:
	./~/clean.sh
