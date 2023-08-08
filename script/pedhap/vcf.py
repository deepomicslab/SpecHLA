"""
Functions for reading VCFs.
"""
import os
import sys
import math
import logging
import itertools
from dataclasses import dataclass
from abc import ABC, abstractmethod
from os import EX_OSFILE, PathLike, set_inheritable
from typing import List, Sequence, Dict, Tuple, Iterable, Optional, Union, TextIO, Iterator

from pysam import VariantFile, VariantHeader, VariantRecord
from read_set import Read, ReadSet
from read_set import ReadSet
from utils import warn_once
logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class VcfError(Exception):
    pass


class VcfNotSortedError(VcfError):
    pass


class PloidyError(VcfError):
    pass


class VcfIndexMissing(VcfError):
    pass


class VcfInvalidChromosome(VcfError):
    pass


class Genotype(object):
    def __init__(self, alleles):
        self.ploidy = 2
        self.alleles = alleles

    def as_vector(self):
        return self.alleles

    def is_homozygous(self):
        return self.alleles[0] == self.alleles[1]


def mendelian_conflict(genotypem, genotypef, genotypec):
    alleles_m = genotypem
    alleles_f = genotypef
    alleles_c = genotypec
    if alleles_c[0] in alleles_m and alleles_c[1] in alleles_f:
        return False
    elif alleles_c[1] in alleles_m and alleles_c[0] in alleles_f:
        return False
    else:
        return True
class VcfVariant:
    """A variant in a VCF file (not to be confused with core.Variant)"""

    __slots__ = ("position", "reference_allele", "alternative_allele")

    def __init__(self, position: int, reference_allele: str, alternative_allele: str):
        """
        Multi-ALT sites are not modelled.
        """
        self.position = position
        self.reference_allele = reference_allele
        self.alternative_allele = alternative_allele

    def __repr__(self):
        return "VcfVariant({}, {!r}, {!r})".format(
            self.position, self.reference_allele, self.alternative_allele
        )

    def __hash__(self):
        return hash((self.position, self.reference_allele, self.alternative_allele))

    def __eq__(self, other):
        return (
            (self.position == other.position)
            and (self.reference_allele == other.reference_allele)
            and (self.alternative_allele == other.alternative_allele)
        )

    def __lt__(self, other):
        return (self.position, self.reference_allele, self.alternative_allele) < (
            other.position,
            other.reference_allele,
            other.alternative_allele,
        )

    def is_snv(self) -> bool:
        return (self.reference_allele != self.alternative_allele) and (
            len(self.reference_allele) == len(self.alternative_allele) == 1
        )

    def normalized(self) -> "VcfVariant":
        """
        Return a normalized version of this variant.

        Common prefixes and/or suffixes between the reference and alternative allele are removed,
        and the position is adjusted as necessary.

        >>> VcfVariant(100, 'GCTGTT', 'GCTAAATT').normalized()
        VcfVariant(103, 'G', 'AAA')
        """
        pos, ref, alt = self.position, self.reference_allele, self.alternative_allele
        while len(ref) >= 1 and len(alt) >= 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]

        while len(ref) >= 1 and len(alt) >= 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos += 1

        return VcfVariant(pos, ref, alt)


@dataclass
class VariantCallPhase:
    block_id: int  # numeric id of the phased block
    phase: List[int]  # alleles representing the phasing. (1, 0) is 1|0
    quality: Optional[int]
    position: int

    def is_homo(self):
        return self.phase[0] == self.phase[1]


class VariantTable:
    """
    For a single chromosome, store variants and their genotypes.
    Each row of this table contains a variant, each column
    contains the genotypes of a single sample.

    chromosome -- chromosome name
    samples -- list of sample names
    """

    def __init__(self, chromosome: str, samples: List[str]):
        self.chromosome = chromosome
        self.samples = samples
        self.phases: List[List[Optional[VariantCallPhase]]] = [[]
                                                               for _ in samples]
        # does this sample phased in prev iteration                                                            
        self.phase_tags: List[bool] = [True for _ in samples]
        self._sample_to_index = {
            sample: index for index, sample in enumerate(samples)}
        self.mendel_cs = []

    def __len__(self) -> int:
        return len(self.phases[0])

    def check_mendel_conflict(self, c, f, m):

        try:
            c_index = self._sample_to_index[c]
            f_index = self._sample_to_index[f]
            m_index = self._sample_to_index[m]
            c_phases = self.phases[c_index]
            f_phases = self.phases[f_index]
            m_phases = self.phases[m_index]
        except KeyError:
                return
        # mendelian_conflicts = []
        for index, (gt_mother, gt_father, gt_child) in enumerate(
                zip(m_phases, f_phases, c_phases)
        ):
            # if (not gt_mother.phase) and (not gt_father.is_none()) and (not gt_child.is_none()):
            if mendelian_conflict(gt_mother.phase, gt_father.phase, gt_child.phase):
                self.mendel_cs.append(gt_child.position)
        logger.info(f"{len(self.mendel_cs)} for contig {self.chromosome}")
        # return mendelian_conflicts

    def write(self, s, out):
        f = open(out, "w")
        s_index = self._sample_to_index[s]
        for p in self.phases[s_index]:
            if p.is_homo():
                # f.write(f"{p.position}\t{p.phase[0]}|{p.phase[1]}\n")
                pass
            elif p.block_id != 0:
                f.write(f"{p.position+1}\t{p.phase[0]}|{p.phase[1]}\t{p.block_id}\n")
            else:
                f.write(f"{p.position+1}\tunphased\n")
        f.close()

    def add_variant(
        self,
        phases: Sequence[Optional[VariantCallPhase]],
    ) -> None:
        if len(phases) != len(self.phases):
            raise ValueError("Expecting as many phases as there are samples")
        for i, phase in enumerate(phases):
            self.phases[i].append(phase)

    def phases_of(self, sample: str) -> List[Optional[VariantCallPhase]]:
        """Retrieve phases by sample name"""
        return self.phases[self._sample_to_index[sample]]

    def num_of_blocks_of(self, sample: str) -> int:
        """ Retrieve the number of blocks of the sample"""
        return len(
            set(
                i.block_id for i in self.phases[self._sample_to_index[sample]] if i is not None)
        )

    def id_of(self, sample: str) -> int:
        """Return a unique int id of a sample given by name"""
        return self._sample_to_index[sample]

    def remove_rows_by_index(self, indices: Iterable[int]) -> None:
        """Remove variants given by their index in the variant list"""
        for i in sorted(indices, reverse=True):
            del self.variants[i]
            for ph in self.phases:
                del ph[i]
        for ph in self.phases:
            assert len(self.variants) == len(ph)
        assert (
            len(self.samples)
            == len(self.phases)
        )

    def subset_rows_by_position(self, positions: Iterable[int]) -> None:
        """Keep only rows given in positions, discard the rest"""
        positions = frozenset(positions)
        to_discard = [i for i, v in enumerate(
            self.variants) if v.position not in positions]
        self.remove_rows_by_index(to_discard)

        # reads_set = Dict[int, Read]
        single_unphased_snps = []
        for variant, genotype, phase in zip(
            self.variants,  self.phases[sample1_index]
        ):
            if genotype.is_homozygous():
                continue
            # if unphased, must not covered by any reads
            if phase is None:
                continue
            if len(genotype.as_vector()) > 2:
                # only use diploid variants
                continue
            if phase is None:
                continue
            if genotype.is_homozygous():
                continue

    # TODO: extend this to polyploid case

    def extend_by_readset(self, s1: str, read_set: ReadSet, side = -1):
        # read_set.finalize()
        try:
            sample1_index = self._sample_to_index[s1]
            sample1_phases = self.phases[sample1_index]
        except KeyError:
            return
        if len(read_set.uncertain_blocks) != 0:
            logger.info(f"{s1}, {len(read_set.uncertain_blocks)} blocks covered but uncertain in this round")
        # if not read_set.contains_phasing_info():
        #     self.phase_tags[sample1_index] = False
        #     logger.info(f"No certain info provided, {s1} unphased in this round")
        #     return
        # phase_info = {}
        finalize_new_block_ids = {}
        ensure_block =[]
        for i in sample1_phases:
            if i.block_id == 0:
                continue
            # if i.block_id in phase_info.keys():
            #     new_block_id = phase_info[i.block_id][0]
            #     need_r = phase_info[i.block_id][1]
            # else:
            new_block_id, need_r = read_set.get_phase_id(i.block_id)
            if new_block_id == 0:
                if i.block_id in finalize_new_block_ids.keys():
                    finalize_new_block_id = finalize_new_block_ids[i.block_id]
                else:
                    finalize_new_block_id = len(finalize_new_block_ids) + 1
                    finalize_new_block_ids[i.block_id] = finalize_new_block_id
            else:
                # phase_info[i.block_id] = [new_block_id, need_r]
                if new_block_id in finalize_new_block_ids.keys():
                    finalize_new_block_id = finalize_new_block_ids[new_block_id]
                else:
                    finalize_new_block_id = len(finalize_new_block_ids) + 1
                    finalize_new_block_ids[new_block_id] = finalize_new_block_id
            i.block_id = finalize_new_block_id
            if need_r:
                t = i.phase[0]
                i.phase[0] = i.phase[1]
                i.phase[1] = t
        for k,v in finalize_new_block_ids.items():
            if k != v:
                ensure_block.append(v)
                ensure_block.append(side)
                break
        return ensure_block



    def phase_with_homo(
        self,
        sample1: str,
        sample2: str,
        prev_ensure_block = None,
        hete_confilic_poses = None,
        default_quality: int = 20,
        mapq: int = 100,
        side: int = 0,
        threshold1: float = 0.1,
        threshold2: float = 0,
    ):
        try:
            sample1_index = self._sample_to_index[sample1]
            sample2_index = self._sample_to_index[sample2]
            sample1_phases = self.phases[sample1_index]
        except KeyError:
            return
        s1_prev_phasing_state = self.phase_tags[sample1_index]
        s2_prev_phasing_state = self.phase_tags[sample2_index]
        # if s1 and s2 unphased in prev round, skip
        if not s1_prev_phasing_state and not s2_prev_phasing_state:
            logging.info(f"Skip phasing {sample1} with {sample2} due to both of them unchanged in prev round")
            # return
        homo_read_set = ReadSet()
        r = Read(mapq, -10101010, threshold1=threshold1, threshold2=threshold2)
        unphase_poses = {}
        for i, phase2 in enumerate(self.phases[sample2_index]):
            phase1 = sample1_phases[i]
            # phase3 = self.phases[3][i]
            if phase1.is_homo():
                continue
            if phase1.position in self.mendel_cs:
                continue
            if phase2.is_homo():
                # if phase1.position == 2402498:
                #     print("xxx")
                if phase1.block_id != 0:
                    o_side = side
                    value = 1
                    if phase2.phase[0] == phase1.phase[1]:
                        o_side = abs(side - 1)
                    if phase2.phase[0] != 0:
                        value = 1
                    r.set_covered_block(
                        phase1.block_id, o_side,phase1.position,value)
                else:
                    # pass
                    # if sample2_index == 2 and phase3.is_homo() and phase2.phase[0] == phase3.phase
                    phase1.block_id = -10101010
                    n_flip = 0
                    if side == 0:
                        if phase1.phase[0] != phase2.phase[0]:
                            # n_flip = 1
                            t = phase1.phase[0]
                            phase1.phase[0] = phase1.phase[1]
                            phase1.phase[1] = t
                    else:
                        if phase1.phase[0] == phase2.phase[0]:
                            # n_flip = 1
                    # unphase_poses[phase1.position] = n_flip
                            t = phase1.phase[0]
                            phase1.phase[0] = phase1.phase[1]
                            phase1.phase[1] = t
        homo_read_set.add_read(r,prev_ensure_block)
        print(homo_read_set.reverse_info,"44444")
        ensure_block = self.extend_by_readset(sample1, homo_read_set, side=side)
        print(ensure_block)
        return homo_read_set.confilict_poses, ensure_block

    def phase_with_hete(
        self,
        sample1: str,
        sample2: str,
        default_quality: int = 20,
        mapq: int = 100,
        side: int = 0,
        threshold1: float = 0.1,
        threshold2: float = 0,
    ):
        try:
            sample1_index = self._sample_to_index[sample1]
            sample2_index = self._sample_to_index[sample2]
            sample1_phases = self.phases[sample1_index]
        except KeyError:
            return
        s1_prev_phasing_state = self.phase_tags[sample1_index]
        s2_prev_phasing_state = self.phase_tags[sample2_index]
        # if s1 and s2 unphased in prev round, skip
        if not s1_prev_phasing_state and not s2_prev_phasing_state:
            logging.info(f"Skip phasing {sample1} with {sample2} due to both of them unchanged in prev round")
            return
        # input_variant_set = set(input_variants)
        heter_read_map: Dict[int, Read] = {}  # maps block_id core.Read objects
        unphase_poses = {}
        for i, phase2 in enumerate(self.phases[sample2_index]):
            # target homo skip,
            phase1 = sample1_phases[i]
            # if phase1.position == 2889794:
            #     print("xxx")
            if phase1.is_homo():
                continue
            if phase2.is_homo():
                continue
            if phase2.block_id == 0:
                continue
            if phase2.quality is None:
                quality = default_quality
            else:
                quality = phase2.quality
            if phase1.position in self.mendel_cs:
                continue
            if phase2.block_id not in heter_read_map:
                r = Read(mapq, -phase2.block_id, threshold1=threshold1, threshold2=threshold2)
                heter_read_map[phase2.block_id] = r
            if phase1.block_id != 0:
                o_side = side
                if {phase1.phase[0], phase1.phase[1]} != {phase2.phase[1], phase2.phase[0]}:
                    if phase2.phase[0] == phase1.phase[1] or phase2.phase[1] == phase1.phase[0]:
                        o_side = abs(side-1)
                    elif phase2.phase[0] == phase1.phase[0] or phase2.phase[1] == phase1.phase[1]:
                        o_side = abs(side)
                else:
                    if phase2.phase[0] == phase1.phase[1]:
                        o_side = abs(side-1)
                heter_read_map[phase2.block_id].set_covered_block(
                    phase1.block_id, o_side, phase1.position)
            else:
                pass
                # if phase1.phase[0] == phase2.phase[1]:
                #     unphase_poses[phase1.position] = 1
                # else:
                #     unphase_poses[phase1.position] = 0
                if {phase1.phase[0], phase1.phase[1]} != {phase2.phase[1], phase2.phase[0]}:
                    continue
                phase1.block_id = -phase2.block_id
                phase1.phase = phase2.phase
        heter_read_set = ReadSet()
        for k, read in heter_read_map.items():
            print(read.covered_blocks, "ccc")
            heter_read_set.add_read(read)
            print(heter_read_set.reverse_info,"ddd")

        self.extend_by_readset(sample1, heter_read_set)
        return heter_read_set.confilict_poses, unphase_poses

    def flip(self, sample):
        try:
            sample_index = self._sample_to_index[sample]
            sample1_phases = self.phases[sample_index]
        except KeyError:
            return
        for i, phase in enumerate(sample1_phases):
            tmp = phase.phase[0]
            phase.phase[0] = phase.phase[1]
            phase.phase[1] = tmp

    def adjust_confilict(self,confilict_poses, sample: str):
        try:
            sample_index = self._sample_to_index[sample]
            sample1_phases = self.phases[sample_index]
        except KeyError:
            return
        for i, phase in enumerate(sample1_phases):
            if phase.position in self.mendel_cs:
                continue
            if phase.position in confilict_poses:
                tmp = phase.phase[0]
                phase.phase[0] = phase.phase[1]
                phase.phase[1] = tmp
            # elif phase.position in unphased_flip_situation[0]:
            #     phase.
class VcfReader:
    """
    Read a VCF file chromosome by chromosome.
    """
    def __init__(
        self,
        path: Union[str, PathLike],
        indels: bool = True,
        phases: bool = False,
        ignore_genotypes: bool = False,
        ploidy: int = None,
    ):
        """
        path -- Path to VCF file
        indels -- Whether to include also insertions and deletions in the list of
            variants.
        ignore_genotypes -- In case of genotyping algorithm, no genotypes may be given in
                                vcf, so ignore all genotypes
        ploidy -- Ploidy of the samples
        """
        # TODO Always include deletions since they can 'overlap' other variants
        self._indels = indels
        self._vcf_reader = VariantFile(os.fspath(path))
        self._path = path
        self._phases = phases
        self._ignore_genotypes = ignore_genotypes
        # intentionally public
        self.samples = list(self._vcf_reader.header.samples)
        self.ploidy = ploidy
        logger.debug("Found %d sample(s) in the VCF file.", len(self.samples))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # follows same structure as for ReadSetReader
        self.close()

    def close(self):
        self._vcf_reader.close()

    @property
    def path(self) -> str:
        return self._vcf_reader.filename.decode()

    def _fetch(self, chromosome: str, start: int = 0, end: Optional[int] = None):
        try:
            records = self._vcf_reader.fetch(chromosome, start=start, stop=end)
        except ValueError as e:
            if "invalid contig" in e.args[0]:
                raise VcfInvalidChromosome(e.args[0]) from None
            elif "fetch requires an index" in e.args[0]:
                raise VcfIndexMissing(
                    "{} is missing an index (.tbi or .csi)".format(self._path)
                ) from None
            else:
                raise
        return records

    def fetch(self, chromosome: str, start: int = 0, end: Optional[int] = None) -> VariantTable:
        """
        Fetch records from a single chromosome, optionally restricted to a single region.

        Return a VariantTable object.
        """
        records = list(self._fetch(chromosome, start=start, end=end))
        return self._process_single_chromosome(chromosome, records)

    def fetch_regions(
        self, chromosome: str, regions: Iterable[Tuple[int, Optional[int]]]
    ) -> VariantTable:
        """
        Fetch records from a single chromosome that overlap the given regions.

        :param regions: a list of start, end tuples (end can be None)
        """
        records = []
        for start, end in regions:
            records.extend(list(self._fetch(chromosome, start=start, end=end)))
        return self._process_single_chromosome(chromosome, records)

    def __iter__(self) -> Iterator[VariantTable]:
        """
        Yield VariantTable objects for each chromosome.

        Multi-ALT sites are skipped.
        """
        for chromosome, records in itertools.groupby(self._vcf_reader, lambda record: record.chrom):
            yield self._process_single_chromosome(chromosome, records)

    @staticmethod
    def _extract_HP_phase(call, pos) -> Optional[VariantCallPhase]:
        hp = call.get("HP")
        if hp is None or hp == (".",):
            return None
        fields = [[int(x) for x in s.split("-")] for s in hp]
        for i in range(len(fields)):
            assert fields[0][0] == fields[i][0]
        block_id = fields[0][0]
        phase = tuple(field[1] - 1 for field in fields)
        return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get("PQ", None))

    @staticmethod
    def _extract_GT_PS_phase(call, pos) -> Optional[VariantCallPhase]:
        if call.get("PS", 0) is None:
            block_id = 0
        else:
            block_id = call.get("PS", 0)
        phase = list(call["GT"])
        return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get("PQ", None), position=pos)

    def _process_single_chromosome(self, chromosome: str, records) -> VariantTable:
        phase_detected = None
        n_snvs = 0
        n_other = 0
        n_multi = 0
        table = VariantTable(chromosome, self.samples)
        prev_position = None
        for record in records:
            if not record.alts:
                continue
            # if len(record.alts) > 1:
            #     # Multi-ALT sites are not supported, yet
            #     n_multi += 1
            #     continue

            pos, ref, alt = record.start, str(record.ref), str(record.alts[0])
            if len(ref) == len(alt) == 1:
                n_snvs += 1
            else:
                n_other += 1
                # if not self._indels:
                #     continue

            if (prev_position is not None) and (prev_position > pos):
                raise VcfNotSortedError(
                    "VCF not ordered: {}:{} appears before {}:{}".format(
                        chromosome, prev_position + 1, chromosome, pos + 1
                    )
                )

            if prev_position == pos:
                warn_once(
                    logger, "Skipping duplicated position %s on chromosome %r", pos + 1, chromosome
                )
                # continue
            prev_position = pos

            # Read phasing information (allow GT/PS or HP phase information, but not both),
            # if requested
            phases = []
            for sample_name, call in record.samples.items():
                phase = None
                for extract_phase, phase_name in [
                    (self._extract_HP_phase, "HP"),
                    (self._extract_GT_PS_phase, "GT_PS"),
                ]:
                    p = extract_phase(call, pos)
                    if p is not None:
                        if phase_detected is None:
                            phase_detected = phase_name
                        elif phase_detected != phase_name:
                            raise MixedPhasingError(
                                "Mixed phasing information in input VCF (e.g. mixing PS "
                                "and HP fields)"
                            )
                        phase = p
                        # check for ploidy consistency and limits
                        phase_ploidy = len(p.phase)

                        if p is None or p.block_id is None or p.phase is None:
                            pass
                        elif self.ploidy is None:
                            self.ploidy = phase_ploidy
                        elif phase_ploidy != self.ploidy:
                            print("phase= {}".format(phase))
                            raise PloidyError(
                                "Phasing information contains inconsistent ploidy ({} and "
                                "{})".format(self.ploidy, phase_ploidy)
                            )
                phases.append(phase)
            table.add_variant(phases)

        logger.debug(
            "Parsed %s SNVs and %s non-SNVs. Also skipped %s multi-ALTs.", n_snvs, n_other, n_multi
        )
        return table


def genotype_code(gt: Optional[Tuple[Optional[int], ...]]) -> Genotype:
    """Return genotype encoded as PyVCF-compatible number"""
    if gt is None:
        result = Genotype([])
    elif any(allele is None for allele in gt):
        result = Genotype([])
    else:
        result = Genotype([allele for allele in gt])  # type: ignore
    return result