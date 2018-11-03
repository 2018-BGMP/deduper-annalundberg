#!/usr/bin/env python

'''deduper pseudocode for BGMP Bi624
Given a SAM file of uniquely mapped reads: remove all PCR duplicates, leaving
only a single copy of each read.
Algorithm should avoid loading everything into memory.
Algorithm should be designed for single-end data with 96 UMIs
Records with matching UMI, chromosome, position and strand are duplicates'''

#intended to be used with samfiles that have been sorted by position

import argparse
import re

def get_arguments():
    parser = argparse.ArgumentParser(description="reads in a SAM file designated as paired_end or not")
    parser.add_argument("-f", "--filename", help="name of file",\
                        required=True, type=str)
    parser.add_argument("-p", "--paired_end", help="optional True (paired end) or False (single end), default is False",\
                        required=False, type=str)
    parser.add_argument("-u", "--umi_file", help="file containing known UMIs",\
                        required=False, type=str)
    return parser.parse_args()

def mk_umi_dict(u_file):
    '''(file) -> dict
    a fxn that will make an UMI dictionary containing only the known valid
    UMIs as keys and value an empty string'''
    umi_dict={}
    with open(u_file) as umis:
        for umi in umis:
            umi=umi.strip('\n')
            umi_dict[umi]=''
    return umi_dict

def check_flag(flag):
    '''(int) -> str
    Fxn takes bitwise flag as int, returns '+' or '-' string to indicate which
    strand read belongs to'''
    if (int(flag) & 16) == 16: #0 in position, (-)strand
        return '-'
    return '+'

def pos_correct(pos, cigar, strand):
    '''(int,str, str)
    Fxn takes in position, cigar string, and strand from a SAM read and adjusts
    the position based on soft clipping and 'N' percieved errors'''
    if strand == '+':
        if 'S' and 'M' in cigar:
            block_cig=re.findall(r'(\d*\D+)', cigar) #isolate cigar pairs as items in a list
            left_s = re.findall(r'(\d+)S', block_cig[0]) #if 1st component is soft clipping grab it
            if len(left_s) == 1:
                pos=pos-int(left_s[0])
    elif strand == '-':
        add_cig=re.findall(r'(\d+)[MND]', cigar) #isolate cigar pairs as items in a list
        soft_cig=re.findall(r'(\d+)S$', cigar) # match all soft clipping
        adj=0
        for num in add_cig:
            adj+=int(num)
        pos=pos+adj-1
        if len(soft_cig)==1:
            pos+=int(soft_cig[0])
    return pos

def check_dup(chrom_dict, umi, pos, strand):
    '''(list, dict)-> bool,dict
    This fxn examines read characteristics and checks if the read is a duplicate
    as compared to chromosome specific dictonary holding written entries.
    Returns a True value for duplicate or False for new entry.'''
    if umi in chrom_dict:
        if strand == chrom_dict[umi][1]:
            if pos == chrom_dict[umi][0]:
                return True
            else:
                return False
        else:
            return False
    else:
        return False

def filter_sam(file,u_flag, u_dict):
    '''(file,bool,dict) -> string
    This fxn iterates through a samfile that is sorted by read position to
    filter out pcr duplicate reads. If u_flag indicates known UMI_dict fxn
    checks if UMI is known before checking for duplicates, if invalid UMI the
    read is discared. uses check_dup to check for duplicate reads, and writes
    only the first entry of each identity to a new_file.sam'''
    chrom=''
    newfile=file.split('/')
    newfile = 'deduped_'+newfile[-1]
    with open(file) as sam, open(newfile, 'w') as new_sam:
        for line in sam:
            if line.startswith('@')==False:
                read_l=line.split('\t')
                if chrom != read_l[2]: #check stored chrom v rname of current read
                    chrom = read_l[2]#set chromosome
                    chrom_dict={} #init chrom-specific dict for kept reads
                qname = read_l[0].split(':')
                umi = qname[len(qname)-1]
                pos, cigar = int(read_l[3]), read_l[5]
                strand = check_flag(read_l[1])
                pos = pos_correct(pos,cigar,strand)
                if u_flag:
                    if umi in u_dict:
                        dup = check_dup(chrom_dict, umi, pos, strand)
                    else:
                        dup = True
                else:
                    dup = check_dup(chrom_dict, umi, pos, strand) #check duplicate and update read_l
                if dup == False:
                    new_sam.write(line) #non duplicate, add entry to new_file
                    chrom_dict[umi]=[pos,strand]
    return newfile


def main():
    '''() -> string
    Uses argparse to get arguments specified in fxn call and runs high level
    functions to return a file name containing non-duplicate sam-reads'''
    args = get_arguments()
    if args.paired_end==True:
        return 'paired end reads not currently supported'
    if args.umi_file is None:
        umi_dict={}
        flag=False
    else:
        umi_dict = mk_umi_dict(args.umi_file)
        flag=True
    newfile_name = filter_sam(args.filename,flag,umi_dict)
    return print(newfile_name)

main()
