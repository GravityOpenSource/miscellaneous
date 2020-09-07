#!/usr/bin/env python3

"""

Alignment2Redis

Author: Xin Wu

GravityOpenSource 2020

Usage: 
    alignment2redis.py  [-s sample_id] [-c config_file]  FASTQ  

Given a sample ID and its FASTQ file, a yaml config file, it uses minimap2 to align reads to a reference genome,
then assign count of reads to corresponding bins of the reference genome, then save and update the redis records.      
The redis records should contain the bins and its depth (coarse resolution for viz) of a genome. If there is more
FASTQ files being parsed, it will update and accummulate the count rather than overwrite it with the same sample id.


Arguments:
    FASTQ        input file in FASTQ format
    
Options:
    -h --help                         Print help
    -s --sample_id=sample_id          Sample ID
    -c --config_file=config_file      Config file [default: alignment2redis.yml] 
    -v --version                      Version
   
""" 

from docopt import docopt
from Bio import SeqIO
import mappy as mp
import numpy as np
import logging, yaml, gzip, math



def load_config(conf_file: str):
    logger.info('loading config file ...')
    with open(conf_file, 'rb') as f:
        conf = yaml.load(f)
        logger.info(conf)
        return conf

# minimap2 pairwise alignment format (PAF):
# query_name, query_len, query_start, query_end, query_target_same_strand, target_name, target_len, target_start, target_end, num_matched, num_base, mapq 
# while mappy alignment object is like:
# q_st  q_en  strand  ctg  ctg_len  r_st  r_en  mlen  blen  mapq  cg:Z:cigar_str
# e.g: 
# 110	902	+	NC_016845.1	5333942	2935103	2935914	750	822	60	tp:A:P	ts:A:.	
#	cg:Z:28M4D1M1D4M1D32M2D1M1D94M1D7M1D15M1D40M1I38M2I12M1D22M1D29M2I94M1D26M1I25M1I14M2I17M1D40M2D71M1D1M2I30M1D25M1D7M2D4M1D10M2D56M1D10M3D28M
# 116	626	-	NC_016845.1	5333942	818646	819163	502	519	60	tp:A:P	ts:A:.	cg:Z:144M1I16M2D1M2D84M1D88M1D67M1I64M2D18M1D26M 
# no matter +/- strand, target end > target start

# return the reference genome length, only regard the first record in FASTA format with gzip compression or not
def get_genome_length(genome_fasta: str):
    genome_length = 0
    if genome_fasta.endswith('.gz'):
        with gzip.open(genome_fasta, 'rt') as handle:
            first_record = True
            for record in SeqIO.parse(handle, 'fasta'): 
                if first_record:
                    genome_length = len(record.seq)
                    first_record = False
                else:
                    break
    else:
        first_record = True
        for record in SeqIO.parse(genome_fasta, 'fasta'):
            if first_record:
                genome_length = len(record.seq)
                first_record = False
            else:
                break
    return genome_length
             

def map(fastq_file: str, conf: dict):
    logger.info('loading config ...')
    genome_name = conf['reference']['genome_name']
    genome_fasta = conf['reference']['fasta_file']
    preset = conf['minimap2']['preset']
    best_n = conf['minimap2']['best_n']
    min_mapq = int(conf['minimap2']['min_mapq'])
    bin_num = int(conf['viz']['bin_num'])

    # reference genome length
    genome_size=get_genome_length(genome_fasta)
    # num of base pairs for a bin
    bin_size = math.ceil(genome_size/bin_num)
    logger.info('reference genome size: {gsize}, bin size: {bsize}, bin num: {bnum}'.format(gsize=genome_size, bsize=bin_size, bnum=bin_num))
    # initiate bins
    bins = np.zeros(bin_num)

    # mapping
    logger.info('starting minimap2 ...')
    aligner = mp.Aligner(genome_fasta, preset=preset, best_n=best_n)
    for name, seq, qual in mp.fastx_read(fastq_file):
        for hit in aligner.map(seq):
            # work around for best_n by only keeping primary mapping, best_n=1 will report 2 hits
            if (hit.is_primary == True) & (hit.mapq >= min_mapq):
                bin_st = math.floor(hit.r_st/bin_size)
                bin_en = math.floor(hit.r_en/bin_size)
                logger.info('--- bins ---')
                print(bin_st, bin_en)
                # for each hit, add 1 for all inbetween bins from hit start bin to end bin (+1 for end not included)
                bins[bin_st : bin_en + 1] += 1
    print(bins)



def main(sample_id: str, fastq_file: str, conf_file: str):
    logger.info(conf_file)
    conf = load_config(conf_file) 
    map(fastq_file, conf) 

    
    


if __name__ == '__main__':
    args = docopt(__doc__, version='alignment2redis version 0.1')
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('a2r')
    logger.info(args)
   
    fastq_file = args['FASTQ']
    sample_id = args['--sample_id']
    conf_file = args['--config_file']
    main(sample_id, fastq_file, conf_file)
