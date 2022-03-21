# commenter/decommenter la bonne ligne pour activer/desactiver le profiler
PROFILER=-pg
#PROFILER=

CFLAGS=-Wall -O2 -c $(PROFILER)

LDFLAGS=-lm $(PROFILER)

all : main test

main : polynome.o main.o
	gcc polynome.o main.o $(LDFLAGS) -o main

test : polynome.o test.o
	gcc polynome.o test.o $(LDFLAGS) -o test

polynome.o : polynome.h polynome.c
	gcc $(CFLAGS) polynome.c

main.o : main.c polynome.h
	gcc $(CFLAGS) main.c

test.o : test.c polynome.h
	gcc $(CFLAGS) test.c

clean : 
	rm -f *.o
	rm -f main
	rm -f test
