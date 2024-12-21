CC=gcc
CFLAGS=-Wall -Wextra -g
LDFLAGS=-lm

all : dos test

dos: dos.c

test: test.c

clean:
	$(RM) dos test
