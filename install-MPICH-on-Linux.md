# install MPICH
```shell
mpich2 
         installing in '/home/mpi/mpich/src'
        procedure:
           ./configure -prefix=/home/mpi/mpich2
           sudo make
           sudo make install
           modify the '~/.bashrc'  
                 +export MPI_ROOT=/home/mpi/mpich2
                 +export PATH=$MPI_ROOT/bin:$PATH
                 +export MANPATH=$MPI_ROOT/man:$MANPATH
           (then you can configure your own code)
           mpicc -o hello hello.c
           mpirun -np 4 ./hello
           
       How to run the hello?
           cd $HOME
           touch .mpd.conf
           chmod 600 .mpd.conf
          (then )
           mpdboot
           cd ~
           touch .mpd.conf
           chmod 600 .mpd.conf
          (then )
           mpirun -np 4 ./hello
 
        (come from internet: http://www.cnblogs.com/liyanwei/archive/2010/04/26/1721142.html  )
```
## sample code for mpich
```c
#include "mpi.h"
#include <stdio.h>
#include <math.h>

int main (int argc, char **argv)
{
int myid, numprocs;
int namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

MPI_Init (&argc, &argv);
MPI_Comm_rank (MPI_COMM_WORLD, &myid);
MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
MPI_Get_processor_name (processor_name, &namelen);
fprintf (stderr, "Hello World! Process %d of %d on %s\n", myid, numprocs, processor_name);
MPI_Finalize ();
return 0;
}
```
## compiled
```shell
$ mpicc -o binary.out main.c -lm -w
```
