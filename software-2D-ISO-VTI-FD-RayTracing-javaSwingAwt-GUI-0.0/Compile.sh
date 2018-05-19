#!/bin/bash

javac MainWindow.java -Xlint:deprecation

javah -d ./lib/ src.seismic.RungeKutta
gcc -shared -fpic -o lib/seismic/libRungeKutta.so \
    -I$JAVA_HOME/include -I$JAVA_HOME/include/linux \
     lib/seismic/src_seismic_RungeKutta.c

java -Djava.library.path=lib/seismic/ MainWindow

echo "Close the Software!"

