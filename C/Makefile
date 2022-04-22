# commenter/decommenter la bonne ligne pour activer/desactiver le profiler
#PROFILER=-pg
PROFILER=

CFLAGS=-Wall -O2 -c $(PROFILER)

LDFLAGS=-lm $(PROFILER)

all : main test

main : polynome.o main.o vect.o
	gcc polynome.o vect.o main.o $(LDFLAGS) -o main

test : polynome.o vect.o test.o
	gcc polynome.o vect.o test.o $(LDFLAGS) -o test

polynome.o : polynome.h polynome.c
	gcc $(CFLAGS) polynome.c

vect.o : vect.h vect.c
	gcc -mavx2 $(CFLAGS) vect.c

main.o : main.c polynome.h
	gcc -mavx2 $(CFLAGS) main.c

test.o : test.c polynome.h
	gcc $(CFLAGS) test.c

clean : 
	rm -f *.o
	rm -f *.ps
	rm -f main
	rm -f test
