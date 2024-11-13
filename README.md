# Scripts, raw results and plotting scripts used for the paper "Signatures of natural selection may indicate a genetic basis for the beneficial effects of oily fish intake in indigenous people from coastal Ecuador"

## data_processing_scripts 
- NMfilter.py script to filter mapped reads for a maximum edit distance 
- beagle2vcf_v3.py script to convert a beagle file with genotype likelihoods to vcf format 

## ohana 

Directory with scripts used to run ohana (job.ohana_structure.sh), plot its results (ohana_plot_q.R) and plot likelihood convergence (plot_ohanalog.R), as well as ohana results files and raw plots.

## pbs 

Directory with output of ANGSD containing Fst and PBS values on windows of 50kb, slid by 10kb (merge2_20M_allchr_snps_ATA_AYM_PEL_fold.w50k.s10k), as well as the scripts used to generate the plots in the paper (plot_angsdPBS.R for genome-wide plots and plot_pbspeaks_AtaAymPEL_w50k.s10k.R for plots zooming into specific genomic regions, with gene annotations). 
