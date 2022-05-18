CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -Wno-unused-function -O2
CPPFLAGS=	-Ihtslib
OBJS=		main.o samview.o vcfview.o bamidx.o bcfidx.o bamshuf.o bam2fq.o tabix.o \
			abreak.o bam2bed.o razf.o razip.o faidx.o bedidx.o pileup.o mapchk.o depth.o genreg.o \
			kthread.o qualbin.o samsort.o peovlp.o bgzip.o
PROG=		htsbox

.SUFFIXES:.c .o
.PHONY:all lib

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

all:$(PROG)

lib:
		cd htslib; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhts.a || exit 1; cd ..

qualbin.o:qualbin.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT $< -o $@

samsort.o:samsort.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT $< -o $@

htsbox:lib $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) -Lhtslib -lhts -lpthread -lz -lm

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG); cd htslib; $(MAKE) clean; cd ..

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE THIS LINE -- make depend depends on it.

abreak.o: boxver.h faidx.h
faidx.o: faidx.h razf.h
fsi.o: htslib/bgzf.h htslib/sam.h htslib/bgzf.h htslib/hts.h
main.o: boxver.h
mapchk.alt.o: faidx.h
mapchk.o: faidx.h
peovlp.o: kvec.h
pileup.o: faidx.h boxver.h
razf.o: razf.h
razip.o: razf.h
