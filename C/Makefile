all : main test

main : polynome.o main.o
	gcc polynome.o main.o -lm -Wall -o main

test : polynome.o test.o
	gcc polynome.o test.o -lm -Wall -o test


polynome.o : polynome.h polynome.c
	gcc -c polynome.c

main.o : main.c polynome.h
	gcc -c main.c

test.o : test.c polynome.h
	gcc -c test.c

clean : 
	rm -f *.o
	rm -f main
	rm -f test