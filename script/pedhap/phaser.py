import sys
from utils import FastaNotIndexedError
from ped_utils import Trio
from typing import Dict, Iterable, List, Optional, TextIO, Tuple
from vcf import VariantCallPhase, VariantTable
from vcf import VcfReader
import logging
from pysam import VariantFile, VariantRecord
logger = logging.getLogger(__name__)

class Phaser(object):
    def __init__(
        self,
        vcf_file: str,
        out_file: TextIO = sys.stdout,
        max_round: int = 5,
        indels: bool = True,
        max_coverage: int = 15,
        tag: str = "PS",
        threshold1: float = 0.1,
        threshold2: float = 0,
    ) -> None:
        super().__init__()
        self.vcf_file = vcf_file
        self.out_file = out_file
        self.vcf_reader = VcfReader(path = self.vcf_file, indels=indels, phases=True)
        self.indels = indels
        self.max_coverage = max_coverage
        self.tag = tag
        self.chromo_variant_table: Dict[str, VariantTable] = {}
        self.chromos = []
        self._variant_file = VariantFile(self.vcf_file)
        self._writer = VariantFile(self.out_file, mode="w", header=self._variant_file.header)
        self._reader_iter = iter(self._variant_file)
        self._unprocessed_record: Optional[VariantRecord] = None
        self._read_vcf()
        self.threshold1 = threshold1
        self.threshold2 = threshold2

    def check_phasing_state(self,chromo):
        v_t = self.chromo_variant_table[chromo]
        for t in v_t.phase_tags:
            if t:
                return True
        return v_t.phase_tags[0]

    def _read_vcf(self):
        logger.info("Reading phased blocks from %r", self.vcf_reader.path)
        for variant_table in self.vcf_reader:
            self.chromo_variant_table[variant_table.chromosome] = variant_table
            self.chromos.append(variant_table.chromosome)

    def phasing_trio_child(self,trio: Trio, chromo):

        child = trio.child
        dad = trio.dad
        mom = trio.mom
        # print(child.id, dad.id, mom.id)

        # self.phasing_duo(child.id, dad.id, chromo, side = 0)
        # self.phasing_duo(child.id, mom.id, chromo, side = 1)
        v_t: VariantTable = self.chromo_variant_table[chromo]
        v_t.check_mendel_conflict(child.id, dad.id, mom.id)
        f_confilict_poses, f_unposes = v_t.phase_with_hete(child.id, dad.id, threshold1=self.threshold1, threshold2=self.threshold2)
        m_confilict_poses, m_unposes = v_t.phase_with_hete(child.id, mom.id, threshold1=self.threshold1, threshold2=self.threshold2)
        fh_confilict_poses, fh_ensure_block = v_t.phase_with_homo(child.id, dad.id,side=0, threshold1=self.threshold1, threshold2=self.threshold2)
        mh_confilict_poses, mh_ensure_block = v_t.phase_with_homo(child.id, mom.id,prev_ensure_block=fh_ensure_block, side=1, threshold1=self.threshold1, threshold2=self.threshold2)


        # f_all = dict(f_unposes.items() + fh_unposes.items())
        # m_all = dict(m_unposes.items() + mh_unposes.items())

        # unphased_flip_situation = {0:[],1:[]}
        # for k,v in f_all.items():
        #     if k in m_all.keys() and v == m_all[k]:
        #         if v == 1:
        #             unphased_flip_situation[1].append(k)
        #         else:
        #             unphased_flip_situation[0].append(k)



        # if f_confilict_poses and m_confilict_poses:
        #     insect_poses = list(set(f_confilict_poses).intersection(set(m_confilict_poses)))
        #     v_t.adjust_confilict(insect_poses,child.id)

        f_m_insect = list(set(f_confilict_poses).intersection(set(m_confilict_poses)))
        f_mh_insect = list(set(f_confilict_poses).intersection(set(mh_confilict_poses)))
        fh_m_insect = list(set(fh_confilict_poses).intersection(set(m_confilict_poses)))
        fh_mh_insect = list(set(fh_confilict_poses).intersection(set(mh_confilict_poses)))

        insects = []
        insects = insects+ f_m_insect
        insects = insects+ f_mh_insect
        insects = insects+ fh_m_insect
        insects = insects+ fh_mh_insect
        v_t.adjust_confilict(insects, child.id)

        # v_t.flip(child.id)


# if fh_confilict_poses and mh_confilict_poses:
    #         insect_poses = list(set(fh_confilict_poses).intersection(set(mh_confilict_poses)))
    #         v_t.adjust_confilict(insect_poses,child.id)


    def phasing_trio_parent(self,trio: Trio, chromo):
        child = trio.child
        dad = trio.dad
        mom = trio.mom
        # print(child.id, dad.id, mom.id)

        self.phasing_duo(dad.id, child.id, chromo, side = 0)
        self.phasing_duo(mom.id, child.id, chromo, side = 0)

    def phasing_duo(self, s1: str, s2: str, chromo, side: int):
        v_t = self.chromo_variant_table[chromo]
        v_t.phase_with_hete(s1, s2)
        v_t.phase_with_homo(s1,s2, side=side)
    def write_simple(self, s1):
        for chromo, v_t in self.chromo_variant_table.items():
            v_t.write(s1,self.out_file)

    def write(self):
        for chromo, v_t in self.chromo_variant_table.items():
            sample_phases: Dict[str, Dict] = dict()
            sample_flip: Dict[str, Dict] = dict()
            for sample in v_t.samples:
                prev_pos = ""
                sample_phases[sample] = {}
                sample_flip[sample] = {}
                for p in v_t.phases[v_t._sample_to_index[sample]]:
                    pos = str(p.position)
                    if p.block_id != 0:
                        if p.block_id in sample_flip[sample]:
                            sample_flip[sample][p.block_id][0] = sample_flip[sample][p.block_id][0] + p.phase[0]
                            sample_flip[sample][p.block_id][1] = sample_flip[sample][p.block_id][1] + p.phase[1]
                        else:
                            sample_flip[sample][p.block_id] = [0,0]
                            sample_flip[sample][p.block_id][0] = p.phase[0]
                            sample_flip[sample][p.block_id][1] = p.phase[1]
                    if  pos in prev_pos:
                        pos = prev_pos + "D"
                        # negative for dup pos
                        sample_phases[sample][pos] = p
                    else:
                        sample_phases[sample][pos] = p
                    # print(prev_pos)
                    prev_pos = pos
            prev_pos = ""
            # print(sample_flip)
            for record in self._record_modifier(chromo):
                pos = str(record.start)
                if not record.alts:
                    continue
                # if len(record.alts) > 1:
                #     # we do not phase multiallelic sites currently
                #     continue
                if pos in prev_pos:
                    # pass
                    pos = prev_pos + "D"
                    # duplicate position, skip it
                    # continue

                for sample in v_t.samples:
                    phase_info = sample_phases[sample][pos]
                    call = record.samples[sample]
                    flip_info = False
                    if call.get("PS", 0) is not None and phase_info.block_id != 0:
                        block_id = phase_info.block_id
                        if sample_flip[sample][block_id][1] < sample_flip[sample][block_id][0]:
                            flip_info = True
                    self._set_PS(
                        call, phase_info, flip_info)
                prev_pos = pos

    def _record_modifier(self, chromosome: str):
        for record in self._iterrecords(chromosome):
            yield record
            self._writer.write(record)

    def _set_PS(
        self,
        call,
        phase: VariantCallPhase,
        flip_info,
    ):
        # assert all(allele in [0, 1] for allele in phase.phase)
        call["PS"] = phase.block_id
        tmp = []
        tmp.append(phase.phase[0])
        tmp.append(phase.phase[1])
        if flip_info:
            tmp[0] = phase.phase[1]
            tmp[1] = phase.phase[0]
            # phase.phase[0] = phase.phase[1]
            # phase.phase[1] = tmp
        call["GT"] = tuple(tmp)
        if phase.is_homo():
            call.phased = False
        if phase.block_id != 0:
            call.phased = True

    def _iterrecords(self, chromosome: str) -> Iterable[VariantRecord]:
        """Yield all records for the target chromosome"""
        n = 0
        if self._unprocessed_record is not None:
            assert self._unprocessed_record.chrom == chromosome
            yield self._unprocessed_record
            n += 1
        for record in self._reader_iter:
            n += 1
            if record.chrom != chromosome:
                # save it for later
                self._unprocessed_record = record
                assert n != 1
                return
            yield record