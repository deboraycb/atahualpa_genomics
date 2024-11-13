#!/bin/bash
#SBATCH --job-name=20Mohana
#SBATCH --output=%A_ohana20M.out
#
# Partition:
#SBATCH --partition=savio
#savio has 164 nodes, 20cores/node, 64GB/node, 0.75 SU/h
#savio2 has 144 (+35) nodes, 24 (or 28) cores/node, 64GB/node, 1.00 SU/h
#savio-bigmem has 4 nodes, 20cores/node, 512GB/node, 1.67 SU/h
#savio2-bigmem has 28nodes, 24 cores/node, 128GB/node, 1.20 SU/h
#
#SBATCH --ntasks=5
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
# Mail type:
#SBATCH --mail-type=end
#
# Mail user:
#SBATCH --mail-user=deboraycb@berkeley.edu
#
# Account:
#SBATCH --account=ac_popgen
#
# Specify Faculty Computing Allowance:
#SBATCH --qos=savio_normal
#
## Command(s) to run:

# load modules used by ohana
module load openblas 
module load r
module load zlib/1.2.11
module load python

## my variables:
OHANADIR=/global/home/users/deboraycb/scratch/atahualpa/analysis/ohana/maf05_3 # dir to write ohana files to
mkdir $OHANADIR
mkdir $OHANADIR/ldprune
mkdir $OHANADIR/lgmfiles
mkdir $OHANADIR/structure
Kmin=3 # min value of k
Kmax=8 # maximum value of k

# PREPARE INPUT FOR STRUCTURE
INPREF=merge2_20M_allchr_snps_maf05
BEAGLE=/global/home/users/deboraycb/scratch/atahualpa/data/beaglemerge/${INPREF}.beagle.gz
# prune down to 100k SNPs (there are 3858409 with maf> 0.05 in total)
zcat /global/home/users/deboraycb/scratch/atahualpa/data/beaglemerge/${INPREF}.beagle.gz| awk 'NR == 1 || NR % 38 == 0' |gzip -c > $OHANADIR/ldprune/${INPREF}.beagle.gz
# convert the beagle file to a lgm file
zcat ${OHANADIR}/ldprune/${INPREF}.beagle.gz | convert bgl2lgm > ${OHANADIR}/lgmfiles/${INPREF}.lgm
# run ohana structure (qpas)
nrep=3
for k in $(seq $Kmin $Kmax); do
    echo "running qpas with k=$k mi=1000 ${nrep}X"
    for i in $(seq $nrep); do
    qpas $OHANADIR/lgmfiles/${INPREF}.lgm -k $k -qo $OHANADIR/structure/${INPREF}_k${k}_mi1000_${i}_q.matrix -fo $OHANADIR/structure/${INPREF}_k${k}_mi1000_${i}_f.matrix -mi 1000 &> $OHANADIR/structure/${INPREF}_k${k}_mi1000_${i}.qpas.log &
    done
done
wait
for k in $(seq $Kmin $Kmax); do
    echo "find seed and iter with highest likelihood for k=$k"
    rm tmp
    for file in $OHANADIR/structure/${INPREF}_k${k}_mi1000_*.qpas.log; do
        awk 'NR==1 {seed=$2}; NR==4 {li=$3}; NF==4 && $3>li {li=$3; mi=$1} END {print seed,mi,li}' $file >> tmp
    done
    seed=$(sort -k3 tmp |head -n1 |cut -f 1 -d " ")
    mi=$(sort -k3 tmp |head -n1 |cut -f 2 -d " ")
    bestli=$(sort -k3 tmp |head -n1 |cut -f 3 -d " ")
    echo "running qpas with k=$k mi=${mi} seed=${seed} bestlik=${bestli}"
    qpas $OHANADIR/lgmfiles/${INPREF}.lgm -k $k -qo $OHANADIR/structure/${INPREF}_k${k}_bestmi${bestmi}_q.matrix -fo $OHANADIR/structure/${INPREF}_k${k}_bestmi${mi}_f.matrix -mi $mi --seed $seed &> $OHANADIR/structure/${INPREF}_k${k}_bestmi${bestmi}.qpas.log &
done
wait
echo "generating structure plots with ohana_plot_q.R and plot likelihood convergence"
for q in $OHANADIR/structure/${INPREF}_k*_bestmi*_q.matrix; do
    pref=$(basename $q _q.matrix)
    Rscript ~/scratch/atahualpa/scripts/plot_ohanalog.R -f $OHANADIR/structure/${pref}.qpas.log &
    ~/scratch/atahualpa/scripts/ohana_plot_q.R $q /global/home/users/deboraycb/scratch/atahualpa/data/atahualpamerge2_sampleinfo375.txt $OHANADIR/poporder.txt &
done 
## estimate population covariances with nemeco
#The -fo output file records the allele frequency inference, which can be used to estimate population covariances with nemeco.
for f in $OHANADIR/structure/${INPREF}_k*_bestmi*_f.matrix; do
    pref=$(basename $f _f.matrix)
    nemeco $OHANADIR/lgmfiles/${INPREF}.lgm $f -co $OHANADIR/structure/${pref}_mi250_c.matrix -mi 250 &> $OHANADIR/structure/${pref}_mi250.nemeco.log
## approximate the estimated covariance matrix into a population tree using cov2nwk
    convert cov2nwk $OHANADIR/structure/${pref}_mi250_c.matrix $OHANADIR/structure/${pref}_mi250_tree.nwk
## produce a visualization of the tree using nwk2svg; 
    convert nwk2svg $OHANADIR/structure/${pref}_mi250_tree.nwk $OHANADIR/structure/${pref}_mi250_tree.svg
done
