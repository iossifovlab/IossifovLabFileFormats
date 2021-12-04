#!/usr/bin/env python

import glob
import sys
import argparse
import datetime
from dataclasses import dataclass
from collections import defaultdict
from typing import List, Dict, Tuple

import pysam

# TODO: Run through VCF validator

# TODO: Store the effect and frequency in the INFO field

# TODO: Handle the limit for the number of open files.
#    Solution 1:
#       If the limit is reached, split the files into batches run a pass over
#       the input for each of the batches.
#    Solution 2:
#       Output one command line for each of the batches. The user can then run
#       these in parralel or sequentially or as cluster jobs etc.
#    Solution 3:
#       If the limit is reach, swith to buffering the output for each file and
#       flushing (open for append; wirte; close) when the buffer is full.

# TODO: Impement three filters for genotype:
#       'all' -- output all positions including the ones where all samples
#            have missing (./.) genotypes.
#       'available' -- remove positions where all samples have missing (./.)
#           genotype [THIS IS THE IMPLEMENTED NOW]
#       'non-reference' -- keep only position that have at least one non
#           reference allele
#       NOTE: In the 'available' and 'non-reference' mode, optimize the set of
#           files whole add_position will be called.

# TODO: Add filters by effect and frequency

# TODO: Set up tests with a small dataset.

# TODO: Add filters by gene!

# TODO: setup.py / add to pypy / iossifovlab in anaconda


def get_file_name(suffix):
    files = list(glob.glob(f"{args.in_dir}/*{suffix}"))
    if len(files) != 1:
        raise Exception(f"Can't find the file with suffix {suffix} "
                        f"in the input directory {args.in_dir}")
    return files[0]


def load_genome_info():
    fai_file_name = get_file_name(".fa.fai")
    genome_name = fai_file_name.split("/")[-1][:-len(".fa.fai")]
    chromosomes = {}
    with open(fai_file_name) as file:
        for line in file:
            chrom, chrom_length, _, _, _ = line.strip("\n\r").split("\t")
            chromosomes[chrom] = int(chrom_length)
    return genome_name, chromosomes


@dataclass
class Person:
    familyId: str
    personId: str
    motherId: str
    fatherId: str
    sex: str
    status: str
    role: str


def load_pedigree() -> Tuple[Dict[str, List[Person]], Dict[str, List[Person]]]:
    families = defaultdict(list)
    persons = defaultdict(list)
    ped_file = get_file_name(".ped")
    with open(ped_file) as file:
        file.readline()  # igore the header
        for line in file:
            cs = line.strip("\n\r").split("\t")
            p = Person(*cs)
            families[p.familyId].append(p)
            persons[p.personId].append(p)
    return persons, families


class VCFOutFile:
    def __init__(self, samples, file_name):
        self.file_name = file_name
        self.samples = samples
        self.samples_set = set(samples)
        self.families = {
            p.familyId for pid in self.samples for p in all_persons[pid]}

    def open(self):
        self.file = open(self.file_name, "w")
        self.write_header()

    def write_header(self):
        print(f"##fileformat=VCFv4.2\n"
              f"##fileDate={datetime.date.today().strftime('%Y%m%d')}\n"
              f"##source=mgtvf_2_0_to_vcf\n"
              f"##reference={genome_id}",
              file=self.file)
        for chrom in used_chromosomes:
            print(f"##contig=<ID={chrom},length={all_chromosomes[chrom]}>",
                  file=self.file)
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              file=self.file)
        print('##FORMAT=<ID=AD,Number=3,Type=Integer,Description="Number of '
              'reads for the reference, alternative, and sometimes other '
              'alleles.">',
              file=self.file)
        print(*'#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT'.split(' '),
              *self.samples, sep="\t", file=self.file)

    def add_position(self, chrom, pos, ref, alt, person_genotypes):
        if not (set(person_genotypes) & self.samples_set):
            return
        gens = [person_genotypes.get(person_id, "./.:.,.,.")
                for person_id in self.samples]

        print(*[chrom, pos, ".", ref, alt, ".", ".", ".", "GT:AD"],
              *gens, sep="\t", file=self.file)

    def close(self):
        self.file.close()


def create_files() -> List[VCFOutFile]:

    arg_families = args.families.split(",") if args.families else None
    arg_persons = args.persons.split(",") if args.persons else None
    if arg_families is not None and arg_persons is not None:
        raise ValueError("create_files should not be called with both "
                         "arg_families and arg_persons set.")
    persons = arg_persons
    if persons is None:
        if arg_families is None:
            persons = list(all_persons.keys())
        else:
            persons = [p.personId for fid in arg_families
                       for p in all_families[fid]]
    r: List[VCFOutFile] = []
    if args.mode == "all":
        r.append(VCFOutFile(persons, f"{args.out_dir}/all.vcf"))
    elif args.mode == "by_family":
        families = defaultdict(list)
        if arg_persons is not None:
            for pid in arg_persons:
                for p in all_persons[pid]:
                    families[p.familyId].append(pid)
        else:
            family_ids = arg_families if arg_families is not None else \
                list(all_families.keys())
            for fid in family_ids:
                for p in all_families[fid]:
                    families[fid].append(p.personId)
        for fid, fpersons in families.items():
            r.append(VCFOutFile(
                fpersons, f"{args.out_dir}/family-{fid}.vcf"))
    elif args.mode == "by_person":
        for pid in persons:
            r.append(VCFOutFile([pid], f"{args.out_dir}/person-{pid}.vcf"))
    else:
        raise ValueError(f"Unnown arg_mode {args.mode}")
    return r


parser = argparse.ArgumentParser(
    description="MGTVF2.0 to VCF",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument('in_dir', default='.', nargs="?",
                    help="The directory where the MGTVF2.0 data is stored.")
parser.add_argument('out_dir', default=".", nargs="?",
                    help="The directory where the VCF file will be sotred.")

parser.add_argument('--mode', default="all",
                    help="The structure of the VCF files: "
                    "'all' - one VCF file (all.vcf);"
                    "'by_family' - a VCF file for each family (family-<familyId>.vcf); "
                    "'by_person' - a VCF file for each person (person-<personId>.vcf).")
parser.add_argument('--families', default=None,
                    help="Comma separated list of family ids to expert.")
parser.add_argument('--persons', default=None,
                    help="Comma separated list of person ids to expert. "
                    "If neither --families nor --persons are provided, "
                    "all persons from the pedigree will be exported.")

parser.add_argument('--region', default=None,
                    help="Genomic region in the form <chr>:<start>-<end>.")


args = parser.parse_args(sys.argv[1:])

all_persons, all_families = load_pedigree()
genome_id, all_chromosomes = load_genome_info()

used_chromosomes = all_chromosomes
if args.region:
    chrom, interval = args.region.split(":")
    start, end = map(int, interval.split("-"))
    used_chromosomes = [chrom]

out_files = create_files()

all_vcf_families = {fid for out_file in out_files for fid in out_file.families}

for out_file in out_files:
    out_file.open()

n_out_poss = 0
with pysam.TabixFile(get_file_name("-family.txt.gz")) as in_file:
    if args.region:
        line_gen = in_file.fetch(chrom, start, end)
    else:
        line_gen = in_file.fetch()
    for line in line_gen:
        # print(line)
        chrom, pos, ref, alt, all_genotypes_str = line.strip(
            "\n\r").split("\t")
        family_genotypes_str = all_genotypes_str.split(";")
        family_genotypes = {}
        for family_genotype_str in family_genotypes_str:
            familyId, members_genotype_str = family_genotype_str.split(" ", 1)
            family_genotypes[familyId] = members_genotype_str
        person_genotypes = {}
        for family_id in all_vcf_families & set(family_genotypes.keys()):
            for member_genotype_str in \
                    family_genotypes[family_id].split(" "):
                pid, role, sex, status, gen, cnts = member_genotype_str.split(
                    ":")
                person_genotypes[pid] = f"{gen}:{cnts}"
        if len(person_genotypes) == 0:
            continue
        for out_file in out_files:
            out_file.add_position(chrom, pos, ref, alt, person_genotypes)
        n_out_poss += 1
        #


for out_file in out_files:
    out_file.close()
