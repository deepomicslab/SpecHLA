CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
#DFLAGS=		-D_USE_KURL
DFLAGS=		-D_USE_KNETFILE
LOBJS=		knetfile.o bgzf.o hts.o vcf.o sam.o tbx.o
INCLUDES=
LIBPATH=

.SUFFIXES:.c .o
.PHONY:lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

libhts.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bgzf.o:bgzf.c bgzf.h knetfile.h khash.h kurl.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_MT -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

hts.o:hts.c bgzf.h bgzf.h khash.h kurl.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) hts.c -o $@

knetfile.o:knetfile.h
kurl.o:kurl.h
vcf.o:vcf.h bgzf.h kstring.h khash.h hts.h
sam.o:sam.h bgzf.h kstring.h hts.h
tbx.o:tbx.h bgzf.h kstring.h hts.h

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ *.a *.so *.dylib

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE THIS LINE -- make depend depends on it.

bgzf.o: bgzf.h knetfile.h
hts.o: bgzf.h hts.h kseq.h khash.h knetfile.h ksort.h
knetfile.o: knetfile.h
kurl.o: kurl.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h
tbx.o: tbx.h hts.h bgzf.h khash.h
vcf.o: kstring.h bgzf.h vcf.h hts.h khash.h kseq.h
