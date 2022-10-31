CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function #-fno-inline-functions -fno-inline-functions-called-once
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o rld0.o sys.o diff.o sub.o unpack.o correct.o dfs.o \
			ksw.o seq.o mag.o unitig.o bubble.o sa.o match.o
PROG=		fermi2
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fermi2:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE THIS LINE -- make depend depends on it.

bubble.o: priv.h mag.h kstring.h kvec.h ksw.h khash.h
correct.o: kvec.h khash.h rld0.h kseq.h ksort.h
dfs.o: kvec.h rld0.h
diff.o: rld0.h kvec.h
ksw.o: ksw.h
mag.o: priv.h mag.h kstring.h kvec.h kseq.h khash.h ksort.h
main.o: fermi2.h rld0.h
match.o: fermi2.h rld0.h kvec.h kseq.h
rld0.o: rld0.h
sa.o: fermi2.h rld0.h kvec.h
seq.o: kstring.h kseq.h
sub.o: rld0.h
unitig.o: kvec.h kstring.h rld0.h mag.h priv.h ksort.h
unpack.o: rld0.h kstring.h kseq.h
