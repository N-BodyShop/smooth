#
# Makefile for smooth.
#
CFLAGS	=   -O3
LIBS	=   -lm

default:	smooth clean

clean:
	rm -f *.o

smooth: main.o kd.o smooth.o
	$(CC) $(CFLAGS) -o smooth main.o kd.o smooth.o $(LIBS)

main.o: main.c kd.h smooth.h
	$(CC) $(CFLAGS) -c main.c

kd.o: kd.c kd.h tipsydefs.h
	$(CC) $(CFLAGS) -c kd.c

smooth.o: smooth.c kd.h smooth.h
	$(CC) $(CFLAGS) -c smooth.c



