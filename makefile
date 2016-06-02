CC=mpicc
CFLAGS= -Wall -std=gnu99
DEPS= grid.h mpi_utils.h time_utils.h
OBJ = main.o grid.o mpi_utils.o time_utils.o
LDFlags= -lm 
INCDIR=/home/support/apps/cports/rhel-6.x86_64/gnu/fftw/3.3.4_mpi/include
LDIR=/home/support/apps/cports/rhel-6.x86_64/gnu/fftw/3.3.4_mpi/lib 

diff: $(OBJ) 
		$(CC) -o $@ $^ -I${INCDIR} -L${LDIR} -lm -lfftw3_mpi $(CFLAGS)

main.o: main.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

grid.o: grid.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

mpi_utils.o: mpi_utils.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

time_utils.o: time_utils.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

clean:
	rm -f $(OBJ) diff


