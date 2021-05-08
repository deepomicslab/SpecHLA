PREFIX=/mnt/delta_WS_1/wangmengyao/HLAPro/script/exon/output/hhh/prefix2

EXE_FERMI2=/mnt/delta_WS_1/wangmengyao/HLAPro/script/exon/../../bin/fermikit/fermi.kit/fermi2
EXE_ROPEBWT2=/mnt/delta_WS_1/wangmengyao/HLAPro/script/exon/../../bin/fermikit/fermi.kit/ropebwt2
EXE_BFC=/mnt/delta_WS_1/wangmengyao/HLAPro/script/exon/../../bin/fermikit/fermi.kit/bfc
GENOME_SIZE=3g
K_EC1=33
K_EC2=55
K_UNITIG=40
K_CLEAN=45
K_TRIM=10
K_MERGE=53
N_THREADS=4

INPUT=cat /mnt/delta_WS_1/wangmengyao/HLAPro/script/exon/output/hhh/extract.fa

SHELL:=/bin/bash
export SHELLOPTS:=errexit:pipefail

all:$(PREFIX).mag.gz

$(PREFIX).ec.fq.gz:
	bash -e -o pipefail -c '$(EXE_BFC) -s $(GENOME_SIZE)  -k $(K_EC1) -t $(N_THREADS) <($(INPUT)) <($(INPUT)) 2> $@.log | gzip -1 > $(PREFIX).ec1.fq.gz'; \
	bash -e -o pipefail -c '$(EXE_BFC) -s $(GENOME_SIZE) -Rk $(K_EC2) -t $(N_THREADS) <($(INPUT)) $(PREFIX).ec1.fq.gz 2>> $@.log | gzip -1 > $@'; \
	rm -f $(PREFIX).ec1.fq.gz

$(PREFIX).flt.fq.gz:$(PREFIX).ec.fq.gz
	$(EXE_BFC) -1s $(GENOME_SIZE) -k $(K_TRIM) -t $(N_THREADS) $< 2> $@.log | gzip -1 > $@

$(PREFIX).flt.fmd:$(PREFIX).flt.fq.gz
	$(EXE_ROPEBWT2) -dNCr $< > $@ 2> $@.log

$(PREFIX).pre.gz:$(PREFIX).flt.fmd
	$(EXE_FERMI2) assemble -l $(K_UNITIG) -m $(K_MERGE) -t $(N_THREADS) $< 2> $@.log | gzip -1 > $@

$(PREFIX).mag.gz:$(PREFIX).pre.gz
	$(EXE_FERMI2) simplify -CSo $(K_CLEAN) -m $(K_MERGE) -T $(K_UNITIG) $< 2> $@.log | gzip -1 > $@

