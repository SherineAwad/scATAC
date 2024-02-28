#! /usr/bin/env python
import sys
import argparse
import math


def modify_bed(infile, outfile): 
    for line in open(infile):
        records = (line.strip()).split("\t")
        new_start = int(records[1]) - 50 
        new_end = int(records[2]) + 50
        new_records = records[0] +'\t' +str(new_start) + '\t' + str(new_end) + '\t' + records[3] +'\t' + records[4] +'\t' + records[5]  
        f = open(outfile, "a+")
        print(new_records, file=f)   
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile')
    parser.add_argument('outfile') 
    args = parser.parse_args()
    bedfile = args.bedfile  
    outfile = args.outfile
    modify_bed(bedfile,outfile) 

if __name__ == '__main__':
    main()

   

