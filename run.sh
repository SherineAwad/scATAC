#!/bin/bash
#SBATCH --job-name run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=10000
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=sherinem@umich.edu
#SBATCH --account=thahoang99 

conda activate py37
#Rscript preprocess.R 17weeks_rep3 
#Rscript preprocess.R 32weeks_rep1
#Rscript preprocess.R 32weeks_rep2 


#Rscript filter.R 17weeks_rep3 100000 30000 1000 1000 500 1.2 2 15
#Rscript filter.R 32weeks_rep1 100000 30000 1000 1000 500 1.2 2 15 
#Rscript filter.R 32weeks_rep2 100000 30000 1000 1000 500 1.2 2 15

#Rscript merge.R merged3samples 17weeks_rep3 32weeks_rep1 32weeks_rep2 

#Rscript analyse.R merged3samples
#Rscript annotateClusters.R merged3samples 

#Rscript removeClusters.R merged3samples
#Rscript diffPeaks.R merged3samples
#Rscript plotPeaks.R merged3samples 
#Rscript addMotifs.R merged3samples
#Rscript plot.R merged3samples


#Rscript toBed.R
#Rscript cat.R
#python modifyBed.py diffPeaks.bed modified_diffPeaks.bed
#bedtools  bedtools getfasta -fi ../REFERENCES/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa.fai  -bed modified_diffPeaks.bed -fo diffPeaks.fasta
#python addFasta.py diffPeaks.fasta allDEGs.csv


#python split.py allpeaks.txt Rod
#python split.py allpeaks.txt Cone 
#python split.py allpeaks.txt RGC
#python split.py allpeaks.txt "Muller glia"

#sort -t , -k 4 Rod.csv -r  > rod_sorted.csv 
#sort -t , -k 4 Cone.csv -r  > cone_sorted.csv 
#sort -t , -k 4 RGC.csv -r  > RGC_sorted.csv
#sort -t , -k 4 Mullerglia.csv -r  > Mullerglia_sorted.csv

Rscript plotPeaks.R merged3samples rod
Rscript plotPeaks.R merged3samples cone
Rscript plotPeaks.R merged3samples RGC
Rscript plotPeaks.R merged3samples Mullerglia
