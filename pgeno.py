# -*- coding: utf-8 -*-

import argparse
import csv
import multiprocessing as mp
from pathlib import Path
import sys
import time
from typing import Iterable, List, Tuple, Union

from Bio import SeqIO
from Bio.Seq import Seq
import regex


def bundle_gen(data_files: Iterable[Path], tag_data: dict) -> Iterable[Tuple[str, Seq, dict]]:
    ''' Generator which yields one sample bundled with the tag data at a time.'''
    for data_file in data_files:
        sample_name = str(data_file.name).rstrip(data_file.suffix)
        with open(data_file, 'r') as fin:
            for record in SeqIO.parse(fin, 'fastq'):
                yield (sample_name, str(record.seq), tag_data)

def get_allele(bundle: Tuple[str, str, dict]) -> Union[Tuple[str, Seq, str], Tuple[None, None, None]]:
    ''' Using the compiled regular expressions in the incoming bundle,
    attempts to identify an allele and returns an outgoing bundle containing
    the identified information'''
    sample_name, seq, tag_data = bundle
    for tag_id, ((f1, f2), (r1, r2)) in tag_data.items():

        forward_upstream = f1.search(seq)
        forward_downstream = f2.search(seq)

        if forward_upstream and forward_downstream:
            _, i = forward_upstream.span()
            j, _ = forward_downstream.span()
            if i <= j:
                allele = seq[i:j]
                return (sample_name, tag_id, allele)

        reverse_upstream = r1.search(seq)
        reverse_downstream = r2.search(seq)
        if reverse_upstream and reverse_downstream:
            _, i = reverse_upstream.span()
            j, _ = reverse_downstream.span()
            if i <= j:
                allele = str(Seq(seq[i:j]).reverse_complement())
                return (sample_name, tag_id, allele)

    return (sample_name, None, None)

def get_tags(data_files: Iterable[Path], tag_files: Iterable[Path], tag_len: int, out_dir: Path, tolerance: int) -> Tuple[dict, dict]:
    ''' Creates and returns two dictionaries, one which contains default
    information for samples which is then used in downstream analysis, and
    another which contains compiled regex information for identifying alleles.
    '''
    data_record = {}
    tag_data = {}
    for tag_file in tag_files:
        for tag_id, tag_info in tags_from_fasta(fasta_file=tag_file, tag_len=tag_len, tolerance=tolerance):
            tag_data[tag_id] = tag_info
            for data_file in data_files:
                sample_name = str(data_file.name).rstrip(data_file.suffix)
                base = {
                    tag_id: {'locus': tag_id}
                }
                try:
                    data_record[sample_name].update(base)
                except KeyError:
                    data_record[sample_name] = base
    return (data_record, tag_data)

def record_intermediate(directory: Path, data_record: dict, allele_set: set) -> None:
    ''' Records data in a longform for potential error and sanity checking.'''
    for fout_name, results in data_record.items():
        fieldnames = ['locus']
        fieldnames.extend(list(allele_set))
        with open(str(directory / f'{fout_name}.csv'), 'w', newline='') as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for tag_id, result in results.items():
                writer.writerow(result)

def tags_from_fasta(fasta_file: Path, tag_len: int, tolerance: int) -> Iterable[tuple]:
    ''' Given a specially annotated `fasta_file`, returns a tuple containing
    `tag_len`-length k-mers compiled into regular expressions which may
    additionally support `tolerance` number of replacements to allow for
    sensing of similar k-mers.
    '''
    with open(fasta_file) as fin:
        for record in SeqIO.parse(fin, 'fasta'):
            upstream_tag_end = record.seq.find('[')
            upstream_tag_start = upstream_tag_end - tag_len
            downstream_tag_start = record.seq.find(']') + 1
            downstream_tag_end = downstream_tag_start + tag_len

            upstream_tag = record.seq[upstream_tag_start:upstream_tag_end]
            downstream_tag = record.seq[downstream_tag_start:downstream_tag_end]
            reverse_upstream_tag = upstream_tag.reverse_complement()
            reverse_downstream_tag = downstream_tag.reverse_complement()

            f1 = regex.compile(f'(?b)({str(upstream_tag)})' + '{s<=' + str(tolerance) + '}')
            f2 = regex.compile(f'(?b)({str(downstream_tag)})' + '{s<=' + str(tolerance) + '}')
            r1 = regex.compile(f'(?b)({str(reverse_downstream_tag)})' + '{s<=' + str(tolerance) + '}')
            r2 = regex.compile(f'(?b)({str(reverse_upstream_tag)})' + '{s<=' + str(tolerance) + '}')

            yield (record.id, ((f1, f2), (r1, r2)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', metavar='./path/to/DATA', type=str, help='path to the folder or file which holds the samples to genotype', required=True)
    parser.add_argument('-l', '--locus', metavar='./path/to/locus.fa', type=str, help='path to the locus file.', required=True)
    parser.add_argument('-o', '--out', metavar='./path/to/OUTPUT', type=str, default='./genotype_data', help='path to the folder where the processed information should be written')
    parser.add_argument('-k', '--kmerlength', metavar='N', type=int, default=10, help='The length of k-mer to be created based on the supplied FASTA files.')
    parser.add_argument('-m', '--minreads', metavar='N', type=int, default=40, help='Minimum number of reads that must be discovered to consider genotyping.')
    parser.add_argument('-g', '--threshold', metavar='N', type=float, default=0.15, help='Fractional representation a sensed haplotype must reach to be considered for genotyping.')
    parser.add_argument('-t', '--threads', metavar='N', type=int, default=0, help='Number of threads to use. 0 to let pgeno decide.')
    parser.add_argument('-z', '--fuzzy', metavar='N', type=int, default=0, help='Number of fuzzy errors (substitutions) allowed when matching k-mers. A higher number is more permissive, but computationally slower.')
    parser.add_argument('--serial', action='store_true', help='Specify this flag to calculate genotypes in serial instead of parallel (will be very slow).')
    args = parser.parse_args()

    data_arg = args.data
    tags_arg = args.loci
    out_dir = args.out
    tag_len = args.kmerlength
    min_reads = args.minreads
    genotype_threshold = args.threshold
    workers = args.threads
    tolerance = args.fuzzy
    serial = args.serial

    if not workers:
        workers = int(mp.cpu_count() * 0.66)

    if workers <= 0 or workers > mp.cpu_count():
       workers = mp.cpu_count() - 1 or 1

    data_files = Path(data_arg)
    tags_files = Path(tags_arg)
    out_dir = Path(out_dir)

    if data_files.is_dir():
        data_files = list(data_files.iterdir())
    else:
        data_files = [data_files]

    if tags_files.is_dir():
        tags_files = list(tags_files.iterdir())
    else:
        tags_files = [tags_files]

    out_dir.mkdir(exist_ok=True)
    intermediate_dir = out_dir / 'intermediate_data'
    intermediate_dir.mkdir(exist_ok=True)

    data_record, tag_data = get_tags(data_files, tags_files, tag_len, out_dir, tolerance)

    start_time = time.time()
    print(f'[{time.ctime()}] Genotyping...')

    bundler = bundle_gen(data_files, tag_data)
    allele_set = set()

    if serial:
        print('Working in serial...')
        for bundle in bundler:
            sample_name, tag_id, allele = get_allele(bundle)
            if tag_id != None and allele != None:
                allele_set.add(allele)
                try:
                    data_record[sample_name][tag_id][allele] += 1
                except KeyError:
                    data_record[sample_name][tag_id][allele] = 1

    else:
        print(f'Working in parallel using {workers} threads...')
        pool = mp.Pool(workers)
        for sample_name, tag_id, allele in pool.imap_unordered(get_allele, bundler, chunksize=1000):
            if tag_id != None and allele != None:
                allele_set.add(allele)
                try:
                    data_record[sample_name][tag_id][allele] += 1
                except KeyError:
                    data_record[sample_name][tag_id][allele] = 1
        pool.close()
        pool.join()

    record_intermediate(intermediate_dir, data_record, allele_set)

    genotypes = {}
    for sample, loci in data_record.items():
        genotypes[sample] = {'sample': sample}
        for record in loci.values():
            locus_name = f"{record['locus']}_genotype"
            locus_dist_name = f"{record['locus']}_distributions"
            genotype = []
            genotype_dist = []
            alleles = record.copy()
            del alleles['locus']
            alleles = sorted(alleles.items(), key=lambda x: x[1], reverse=True)
            total = sum(a[1] for a in alleles)
            if total < min_reads:
                continue
            alleles = [a for a in alleles if a[1] / total > genotype_threshold]

            if len(alleles) < 1:
                genotype.append('None')
                genotype_dist.append('')
            elif len(alleles) == 1:
                allele, count = alleles[0]
                if allele == '':
                    allele = '-'
                genotype = [allele, allele]
                genotype_dist.append(str(count))
            elif len(alleles) > 1:
                for allele, count in alleles[:2]:
                    if allele == '':
                        allele = '-'
                    genotype.append(allele)
                    genotype_dist.append(str(count))

            geno_str = '/'.join(genotype)
            geno_dist_str = '/'.join(genotype_dist)
            genotypes[sample][locus_name] = geno_str
            genotypes[sample][locus_dist_name] = geno_dist_str

    with open(str(out_dir / 'genotypes.csv'), 'w', newline='') as fout:
        headers = ['sample']
        loci_geno = [f'{locus}_genotype' for locus in tag_data.keys()]
        loci_dist = [f'{locus}_distributions' for locus in tag_data.keys()]
        for locus, locus_dist in zip(loci_geno, loci_dist):
            headers.extend([locus, locus_dist])
        writer = csv.DictWriter(fout, fieldnames=headers, extrasaction='ignore')
        writer.writeheader()
        for row in genotypes.values():
            writer.writerow(row)

    finish_time = time.time() - start_time
    print(f'[{time.ctime()}] Finished genotyping in {finish_time:.2f}s.')
