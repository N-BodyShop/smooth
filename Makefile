#
# Makefile for smooth.
#
CFLAGS	=   -O3
LIBS	=   -lm

default:	smooth

clean:
	rm -f *.o

smooth: main.c kd.c smooth.c kd.h smooth.h tipsydefs.h
	$(CC) $(CFLAGS) -o smooth main.c kd.c smooth.c $(LIBS)


