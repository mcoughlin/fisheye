CC= cc
CFLAGS= -O -g -Wall
ULOC= .

fisheye: fisheye.o fitlm.o linsolve.o
	$(CC) $(CFLAGS) fisheye.o -o fisheye fitlm.o linsolve.o -lm
