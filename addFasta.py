#! /usr/bin/env python
import sys
import argparse
import math

def read_fasta(fasta_file): 
    seqs= {}
    count = 0 
    with open(fasta_file, 'rb') as f:
        while True:
           line1 = f.readline()
           count +=1 
           if not line1:
               break
           line2 = f.readline() 
           sid = (line1.strip())
           sid = sid.decode("utf-8")
           sid =sid.replace('>','')
           sid =sid.replace("(.)", '')
           seq = line2.strip()
           seq = seq.decode("utf-8")
           ids = sid.split(":")
           chrs = ids[0] 
           pos = ids[1]
           pos = pos.split("-")
           pos1 =pos[0]
           pos2 = pos[1] 
           id =chrs+"-"+pos1+"-"+pos2
           seqs[id]= seq
    return seqs 

def read_csv(infile, seqs): 
    records ={} 
    count = 0
    f = open("allpeaks.txt", "a+")
    for line in open(infile):
        if count == 0: 
              count+=1
              print(line.strip(),",","Sequence", file=f) 
              continue 
        records = (line.strip()).split(",")
        pos= records[0].replace('"', '')
        pos = pos.split("-") 
        new_start = int(pos[1] ) - 50
        new_end = int(pos[2]) + 50
        modified_pos = pos[0]+"-"+str(new_start)+"-"+str(new_end)
        seq = seqs.get(modified_pos)
        print(line.strip(),",",seq.strip(), file=f)   
        
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file')
    parser.add_argument('infile') 
    args = parser.parse_args()
    fasta_file = args.fasta_file  
    infile = args.infile
    seqs ={}
    seqs = read_fasta(fasta_file) 
    read_csv(infile,seqs) 

if __name__ == '__main__':
    main()

   

