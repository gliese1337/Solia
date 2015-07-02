CC = g++
CFLAGS = -Wall -O3 -I inc/

nbody: nbody.o nbodyio.o evolve.o
	${CC} ${CFLAGS} obj/evolve.o obj/nbodyio.o obj/nbody.o -o nbody

nbody.o: src/nbody.cpp inc/nbody.h inc/nbodyio.h inc/evolve.h
	${CC} ${CFLAGS} -c src/nbody.cpp -o obj/nbody.o

nbodyio.o: src/nbodyio.cpp inc/nbody.h inc/nbodyio.h
	${CC} ${CFLAGS} -c src/nbodyio.cpp -o obj/nbodyio.o

evolve.o: src/evolve.cpp inc/nbody.h inc/evolve.h
	${CC} ${CFLAGS} -c src/evolve.cpp -o obj/evolve.o

clean:
	rm obj/*.o