import argparse
from phaser import Phaser
import logging
from typing import List 
from pedigree import (
    mendelian_conflict,
)
import ped
import ped_utils
from ped_utils import Trio
logger = logging.getLogger(__name__)

def find_mendelian_conflicts(trios, variant_table):
    mendelian_conflicts = set()
    for trio in trios:
        genotypes_mother = variant_table.genotypes_of(trio.mother)
        genotypes_father = variant_table.genotypes_of(trio.father)
        genotypes_child = variant_table.genotypes_of(trio.child)

        for index, (gt_mother, gt_father, gt_child) in enumerate(
            zip(genotypes_mother, genotypes_father, genotypes_child)
        ):
            if (not gt_mother.is_none()) and (not gt_father.is_none()) and (not gt_child.is_none()):
                if mendelian_conflict(gt_mother, gt_father, gt_child):
                    mendelian_conflicts.add(index)
    return mendelian_conflicts

# def phasing_trio_child(phaser: Phaser, trio: Trio, chromo: str):
#     phaser.phasing_trio_child(chromo, trio)


def up_to_down(all_trios: List[Trio], phaser: Phaser, chromo):
    # find top level trio and phsing child
    top_level_trios = ped_utils.get_top_level_trios(all_trios)
    # t = top_level_trios[0]
    # phaser.phasing_trio_child(t)
    for t in top_level_trios:
        phaser.phasing_trio_child(t, chromo)
    next_level_trios = ped_utils.get_next_level_trios(
        all_trios, top_level_trios)
    while len(next_level_trios) != 0:
        for t in next_level_trios:
            phaser.phasing_trio_child(t, chromo)
        next_level_trios = ped_utils.get_next_level_trios(
            all_trios, next_level_trios)


def down_to_up(all_trios: List[Trio], phaser: Phaser, chromo):
    # find bottom level trio and phsing child
    bottom_level_trios = ped_utils.get_bottom_level_trios(all_trios)
    # t = top_level_trios[0]
    # phaser.phasing_trio_child(t)
    # phaser.write_phased_result(t.child.id, "/home/caronkey/Documents/cityu/pedhap/test/test2.out")
    for t in bottom_level_trios:
        phaser.phasing_trio_parent(t, chromo)
    prev_level_trios = ped_utils.get_prev_level_trios(
        all_trios, bottom_level_trios)
    while len(prev_level_trios) != 0:
        for t in prev_level_trios:
            phaser.phasing_trio_parent(t, chromo)
        prev_level_trios = ped_utils.get_prev_level_trios(
            all_trios, prev_level_trios)


def iter_phase(all_trios: List[Trio], phaser: Phaser):
    up_to_down(all_trios, phaser)
    # down_to_up(all_trios, phaser)


def main():
    parser = argparse.ArgumentParser("trio phase")
    parser.add_argument(
        '-v', help='merged VCF file', required=True, dest='vcf_file')
    parser.add_argument(
        '-p', help='pedigree file', required=True, dest='ped_file')
    parser.add_argument(
        '-o', help='out phased vcf file', required=True, dest='out_file')
    parser.add_argument(
        '--threshold1', help='merge conflict blocks threshold', required=False,type=float, default=0.1)
    parser.add_argument(
        '--threshold2', help='merge unconflict blocks threshold', required=False, type=float, default=0)
    parser.add_argument(
        '--max_round', help='max phasing iter times, if not given, decided by program', 
        required=False, default=0, type=int,
        dest='max_round')
    parser.add_argument(
        '--support_c', help='max phasing iter times, if not given, decided by program',
        required=False, default=6, type=int,
        dest='support_c')
    args = parser.parse_args()

    phaser = Phaser(vcf_file=args.vcf_file, out_file=args.out_file, max_round = int(args.max_round), threshold1=args.threshold1, threshold2=args.threshold2)
    families = ped.open_ped(args.ped_file)
    for f in families:
        all_trios = ped_utils.get_trios(f)
        for chromo in phaser.chromos:
            pass
            # up_to_down(all_trios, phaser, chromo)
            # down_to_up(all_trios, phaser, chromo)
            # while phaser.check_phasing_state(chromo):
            up_to_down(all_trios, phaser, chromo)
                # down_to_up(all_trios, phaser, chromo)
    phaser.write()
    # phaser.write_simple("s0210-1_FDHG190451805-1a")
        # phaser.write_phased_result("s0210-1_FDHG190451805-1a", "/home/caronkey/Documents/cityu/pedhap/test/test2.out")
    # run_pedhap(variant_file=args.vcf_file,
    #            ped=args.ped_fle, output=args.out_file)


if __name__ == "__main__":
    main()
