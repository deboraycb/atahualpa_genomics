#install.packages("RCurl")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("Gviz")
library(Gviz) 
library(ensembldb)
library(ggplot2)

#========================================================= 
# TOP REGIONS w50k

fname <- "merge2_20M_allchr_snps_ATA_AYM_PEL_fold.w50k.s10k"
tab <- read.table(fname, header = T)

tab <- tab[order(tab$chr),]
tab$midPos <- tab$midPos/1000
tab$plotPos <- tab$midPos
chr <- 1
for (l in 1:nrow(tab)){
    if (tab$chr[l] != chr){
        chr <- tab$chr[l]
        tab$plotPos[tab$chr==chr] <- (tab$plotPos[tab$chr==chr] + tab$plotPos[l-1])
    }
}
tab$chr <- as.factor(tab$chr)

nwin <- nrow(tab)
perc.1 <- sort(tab$PBS0, decreasing=T)[round(nwin/1000)]
tabsortPBS <- tab[order(tab$PBS0, decreasing=T),c(2,3,4,8)]
tabtop <- tabsortPBS[tabsortPBS$PBS0>perc.1,]
regions <- list()
for(chr in 1:22){regions[chr] <- 0}
tabtop$otherwin <- ""
tabtop$notherwin <- 0
tabtop$maxNsites <- 0
for(l in 1:nrow(tabtop)){
    chr <- tabtop$chr[l]
    pos <- round(tabtop$midPos[l])
    if ( ! pos %in% regions[[chr]]){
        replines <- which(tabtop$midPos[tabtop$chr==chr] > pos-500 & tabtop$midPos[tabtop$chr==chr] < pos+500)
        tabtop$notherwin[l] <- length(replines)
        tabtop$maxNsites[l] <- max(tabtop$Nsites[replines])
        tabtop$otherwin[l] <- paste(round(tabtop$midPos[tabtop$chr==chr][replines]), collapse=",")
        regions[[chr]] <- c(regions[[chr]], round(pos-500):round(pos+500))
    }
}
dim(tabtop)
sum(tabtop$notherwin)
head(tabtop[tabtop$notherwin>0,])
tabtop[tabtop$notherwin>0,-5]

nwin <- nrow(tab)
top.1percPBS0 <- sort(tab$PBS0, decreasing=T)[round(nwin/1000)]
#========================================================= 
# PLOT PEAKS

getregion <- function(tab,chr,start,end,pop1,pop2,pop3){
    subset <- tab[tab$chr==chr & tab$midPos > start/1000 & tab$midPos < end/1000,] 
    itrack <- IdeogramTrack(genome="hg19", chromosome=paste0("chr",chr))
    gtrack <- GenomeAxisTrack()
    pbstrack <- DataTrack(data=subset$PBS0,start=(subset$midPos-25)*1000,
                        end=(subset$midPos+25)*1000, chromosome=paste0("chr",chr), genome="hg19", name="PBS", baseline=top.1percPBS0, col.baseline="blue", lty.baseline=2) 
    fst1track <- DataTrack(data=subset$Fst01,start=(subset$midPos-25)*1000,
                        end=(subset$midPos+25)*1000, chromosome=paste0("chr",chr), genome="hg19", name=paste0("Fst ",pop1,"-",pop2), baseline=0.044464, col.baseline="black", col="gray")
    fst2track <- DataTrack(data=subset$Fst02,start=(subset$midPos-25)*1000,
                        end=(subset$midPos+25)*1000, chromosome=paste0("chr",chr), genome="hg19", name=paste0("Fst ",pop1,"-",pop3), col="gray", baseline=0.040006, col.baseline="black")
    fst3track <- DataTrack(data=subset$Fst12,start=(subset$midPos-25)*1000,
                        end=(subset$midPos+25)*1000, chromosome=paste0("chr",chr), genome="hg19", name=paste0("Fst ",pop2,"-",pop3), col="gray", baseline=0.016423, col.baseline="black")
    rsGenes <- UcscTrack(genome = "hg19", chromosome = paste0("chr",chr), 
                        track = "refSeq", table="refGene",
                        from = start, to = end,
                        trackType = "GeneRegionTrack", 
                        rstarts = "exonStarts", rends = "exonEnds",
                        gene = "name2", symbol = "name", 
                        transcript = "name2", strand = "strand",
                        name = "RefSeq Genes")
    return(list(itrack,gtrack,pbstrack,fst1track,fst2track,fst3track,rsGenes))
}

plotmytracks <- function(tracks,from,to,colap=TRUE,limy=TRUE, extendy=1){
    if(limy){
        lims=range(sapply(tracks[4:6], function(x) range(x[1]@data)))
        lims[2] <- lims[2]*extendy
    }
    else lims=NULL
    plotTracks(tracks, from=from, to=to,
               transcriptAnnotation = "gene",
               collapseTranscripts = colap, shape="arrow",
               background.title = "gray60",
               ylim=lims)
}

start=20e6;end=30e6
toplot <- getregion(tab,1,start,end, "Ata", "Aym", "PEL")
png("peakchr1_26Mb.png", height=3.5, width=5, units="in", res=300)
plotmytracks(toplot, from=25.71e6, to=26.13e6, limy=TRUE, colap=TRUE)
dev.off()
start=65e6;end=75e6
toplot <- getregion(tab,1,start,end)
png("peakchr1_70Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=start, to=end, limy=TRUE)
dev.off()
png("peakchr1_70Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=71e6, to=73e6, limy=TRUE)
dev.off()
start=100e6;end=150e6
toplot <- getregion(tab,1,start,end)
png("peakchr1_120Mb_WARS2TBX15.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=119e6, to=121e6, colap=TRUE)
dev.off()
start=150e6;end=180e6
toplot <- getregion(tab,1,start,end, "Ata", "Aym", "PEL")
png("peakchr1_172Mb.png", height=3, width=5, units="in", res=300)
plotmytracks(toplot, from=171.5e6, to=172.5e6)
dev.off()
png("peakchr1_155Mb_all.png", height=3, width=5, units="in", res=300)
plotmytracks(toplot[1:6], from=155e6, to=159e6, colap=TRUE)
dev.off()
png("peakchr1_155Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=155.1e6, to=156e6, extendy=1.1)
dev.off()
png("peakchr1_155Mb_zoom2.png", height=5, width=5, units="in", res=300)
start=158.5e6;end=159e6
plotmytracks(toplot, from=start, to=end)
dev.off()

start=100e6;end=200e6
toplot <- getregion(tab,2,start,end,"Ata", "Aym", "PEL")
png("peakchr2_190Mb_all.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=188e6, to=192e6)
dev.off()
png("peakchr2_190Mb_zoom1.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=188.5e6, to=190.5e6, colap=FALSE)
dev.off()
png("peakchr2_190Mb.png", height=4, width=5, units="in", res=300)
plotmytracks(toplot, from=190.45e6, to=191.15e6, extendy=1.5)
dev.off()
png("peakchr2_133Mb.png", height=3.5, width=5, units="in", res=300)
plotmytracks(toplot, from=132.65e6, to=133.1e6)
dev.off()
png("peakchr2_143Mb.png", height=3.5, width=5, units="in", res=300)
plotmytracks(toplot, from=140.5e6, to=144.5e6, extendy=1.3)
dev.off()
start=1e6;end=50e6
toplot <- getregion(tab,2,start,end,"Ata", "Aym", "PEL")
png("peakchr2_17Mb.png", height=3.5, width=5, units="in", res=300)
plotmytracks(toplot, from=16.35e6, to=17.55e6, extendy=1.3)
dev.off()

start=100e6;end=200e6
toplot <- getregion(tab,6,start,end, "Ata", "Aym", "PEL")
png("peakchr6_Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=130e6, to=140e6)
dev.off()
png("peakchr6_Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=100e6, to=101e6)
dev.off()

start=50e6;end=70e6
toplot <- getregion(tab,7,start,end,"Ata", "Aym", "PEL")
png("peakchr7_57Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=57e6, to=58e6)
dev.off()
png("peakchr7_100Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=100e6, to=101e6)
dev.off()
start=95e6;end=105e6
toplot <- getregion(tab,7,start,end)
png("peakchr7_100Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot[1:6], from=99e6, to=102e6)
dev.off()
png("peakchr7_100Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=100e6, to=101e6)
dev.off()

start=40e6;end=50e6
toplot <- getregion(tab,8,start,end,"Ata", "Aym", "PEL")
png("peakchr8.png", height=1.5, width=5, units="in", res=300)
plotTracks(toplot[1:3], from=start, to=end,
               background.title = "gray60")
dev.off()
png("peakchr8_40Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=35e6, to=45e6)
dev.off()

start=22e6;end=32e6
toplot <- getregion(tab,10,start,end)
png("peakchr10_28Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot[1:6], from=start, to=end, limy=TRUE)
dev.off()
png("peakchr10_28Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=29e6, to=30e6, limy=TRUE)
dev.off()
start=5e6;end=15e6
toplot <- getregion(tab,10,start,end)
png("peakchr10_10Mb.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot[1:6], from=start, to=end, limy=TRUE)
dev.off()
png("peakchr10_10Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=9e6, to=10.5e6, limy=TRUE)
dev.off()
start=100e6;end=110e6
toplot <- getregion(tab,10,start,end,"Ata", "Aym", "PEL")
png("peakchr10_105Mb.png", height=6, width=5, units="in", res=300)
plotmytracks(toplot[1:6], from=102e6, to=108e6, limy=TRUE)
dev.off()
png("peakchr10_105Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=104100000, to=105300000, limy=TRUE)
dev.off()

start=60e6;end=70e6
toplot <- getregion(tab,11,start,end)
png("peakchr11_61Mb_FADS.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=61.4e6, to=61.9e6, colap=TRUE, limy=TRUE)
dev.off()

start=30e6;end=50e6
toplot <- getregion(tab,17,start,end,"Ata", "Aym", "PEL")
png("peakchr17_36Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=35e6, to=37e6, limy=TRUE)
dev.off()
png("peakchr17_45Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=45e6, to=46e6, limy=TRUE)
dev.off()

start=30e6;end=50e6
toplot <- getregion(tab,16,start,end,"Ata", "Aym", "PEL")
png("peakchr16_34Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=33.5e6, to=34.5e6, limy=TRUE)
dev.off()
png("peakchr16_47Mb_zoom.png", height=5, width=5, units="in", res=300)
plotmytracks(toplot, from=46e6, to=47e6, limy=TRUE)
dev.off()

