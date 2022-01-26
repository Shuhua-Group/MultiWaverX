FLAGS=-fopenmp
all: MultiWaverX
	mkdir -p ../bin
	cp MultiWaverX ../bin/

MultiWaverX: MultiWaverX.cpp ancSegs.o argparse.o bootstrap.o exception.o models.o
	g++ MultiWaverX.cpp ancSegs.o argparse.o bootstrap.o exception.o models.o -o MultiWaverX ${FLAGS} -O3

ancSegs.o: ancSegs.cpp ancSegs.h
	g++ -c ancSegs.cpp -O3

argparse.o: argparse.cpp argparse.h
	g++ -c argparse.cpp -O3

bootstrap.o: bootstrap.cpp bootstrap.h
	g++ -c bootstrap.cpp -O3

exception.o: exception.cpp exception.h
	g++ -c exception.cpp -O3

models.o: models.cpp models.h
	g++ -c models.cpp -O3

clean:
	rm -rf *.o MultiWaverX

LDFLAG=-lnlopt -lm
