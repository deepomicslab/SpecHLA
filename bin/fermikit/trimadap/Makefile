CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

all:trimadap-mt

trimadap-mt:trimadap-mt.c ksw.c kthread.c kseq.h ksw.h
		$(CC) $(CFLAGS) -pthread ksw.c kthread.c trimadap-mt.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session* trimadap-mt
