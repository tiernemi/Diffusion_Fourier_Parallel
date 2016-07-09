CC=mpicc
CFLAGS= -Wall -std=gnu99
DEPS= grid.h mpi_utils.h
OBJ = main.o grid.o mpi_utils.o
LDFlags= -lm 
INCDIR=
LDIR=

diff: $(OBJ) 
		$(CC) -o $@ $^ -I${INCDIR} -L${LDIR} -lm -lfftw3 -lfftw3_mpi $(CFLAGS)

main.o: main.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

grid.o: grid.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

mpi_utils.o: mpi_utils.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS) -I${INCDIR} -L${LDIR} 

animation: diff
	rm -f out*
	mpirun -n 2 ./diff -n 100 -m 100 -p -i 1000 -t 2
	gnuplot -e "load 'animation.gnu'"

clean:
	rm -f $(OBJ) diff


