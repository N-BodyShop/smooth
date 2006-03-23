#
# Makefile for smooth.
#
CFLAGS	=   -O
LIBS	=   -lm

default:	smooth

clean:
	rm -f *.o

smooth: main.o kd.o smooth.o
	$(CC) $(CFLAGS) -o smooth $^ $(LIBS)

kd.o: kd.h tipsydefs.h
main.o: kd.h smooth.h
smooth.o: smooth.h kd.h
