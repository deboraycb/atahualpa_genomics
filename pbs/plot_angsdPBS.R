library(ggplot2)
library(scales)

fname <- "merge2_20M_allchr_snps_ATA_AYM_PEL_fold.w50k.s10k"

tab <- read.table(fname, header = T)

tab <- tab[order(tab$chr),]
tab$midPos <- tab$midPos/1000
tab$plotPos <- tab$midPos- tab$midPos[1]
chr <- 1
for (l in 1:nrow(tab)){
    if (tab$chr[l] != chr){
        chr <- tab$chr[l]
        tab$plotPos[tab$chr==chr] <- tab$plotPos[tab$chr==chr] + tab$plotPos[l-1]- tab$midPos[1]
    }
}
chrmids <- rep(0,22)
chrmids[1] <- tail(tab$plotPos[tab$chr == 1], n=1) / 2
for (chr in 2:22){
    chrmids[chr] <- tail(tab$plotPos[tab$chr == chr-1], n=1) + tail(tab$midPos[tab$chr == chr], n=1)/2
}
tab$chr <- as.factor(tab$chr)
nwin <- nrow(tab)

top.1percPBS0 <- sort(tab$PBS0, decreasing=T)[round(nwin/1000)]

chrcol=rep(c("grey25","grey50"),11)
ylim0=range(tab$PBS0)
ggplot(tab)+
    geom_point(aes(plotPos, PBS0, col=chr), alpha=0.5)+
    geom_hline(aes(yintercept=top.1percPBS0), linetype="dashed", colour="blue", guide=FALSE)+
    scale_color_manual(values=chrcol, guide=FALSE)+
    scale_x_continuous(breaks=chrmids, labels=1:22, expand=expansion(mult=c(0.005, 0.005)))+
    ylim(ylim0[1],ylim0[2])+
    xlab("Chromosome")+
    ylab("PBS")+
    theme_classic()
ggsave("merge2_20M_allchr_snps_ATA_AYM_PEL.w50k.s10k_chronX.png", h=1.5, w=9)
