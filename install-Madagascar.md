# Install Madagascar
```shell
Madagascar   (wrong install,so make installing again,see in the following part)
           sudo apt-get install freeglut3-dev g++ gfortran libgd2-xpm-dev libglew1.5-dev
                            there are some problems about installing (suggested packages:
                                                                         libxcb-doc
                                                                            .....
           sudo apt-get install libjpeg62-dev libx11-dev
           sudo apt-get install libxaw7-dev libnetpbm10-dev swig python-dev python-scipy python-numpy
           sudo apt-get install libtiff4-dev scons units libblas-dev
           sudo apt-get install libcairo2-dev libavcodec-dev libplplot-dev
                            there are some problem too.
           cp -a /Desktop/madagascar /usr
           cd madagascar
           scons    (problem?)
           scons install
           vi .bashrc +

---------------------------------------------------------------------------------------
Madagascar(try installing again,follow the 'INSTALL.txt')
           
             (c complier and Python interpreter have not ready yet)
             cd /usr/madagascar
             ./configure --prefix=/usr/madagascar/install 
             ./configure --prefix=/usr/madagascar/install \API fortran-90
             make   ...Done
             make install
             source /usr/madagascar/install/share/madagascar/etc/env.sh
             #source /usr/madagascar/install/share/madagascar/etc/env.csh

             install is Done,the following code is try to test the installing successfully or not.
             (sfin
	      sfattr
	      sfspike
	      sfbandpass
	      sfwiggle)
```
