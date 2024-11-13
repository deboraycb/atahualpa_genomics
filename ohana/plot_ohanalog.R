library("optparse")
library("ggplot2")

option_list = list(make_option(c("-f", "--file"), type="character", default=NULL, 
	       help="dataset file name", metavar="character"))
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

file <- opt$file

log <- read.table(file, skip=2, header=T, as.is=T, fill=T)

for(col in 1:ncol(log)){
	log[,col] <- as.numeric(log[,col])
}

plotit <- function(toplot, outpref="ohanalog"){
	ggplot(toplot)+
		#geom_point(aes(iter, log_likelihood), color="red")+
		geom_line(aes(iter, log_likelihood), color="red")
	ggsave(paste0(outpref,".llik.png"), h=3, w=6)
	ggplot(toplot)+
		#geom_point(aes(iter, delta.lle), color="blue")+
		geom_line(aes(iter, delta.lle), color="blue")
	ggsave(paste0(outpref,".deltal.png"), h=3, w=6)
}

plotit(log[1:500,], outpref=paste0(file, ".000-500"))
plotit(log[10:500,], outpref=paste0(file, ".010-500"))
plotit(log[30:500,], outpref=paste0(file, ".030-500"))
plotit(log[50:500,], outpref=paste0(file, ".050-500"))
plotit(log[70:500,], outpref=paste0(file, ".070-500"))
plotit(log[90:500,], outpref=paste0(file, ".090-500"))
plotit(log[110:500,], outpref=paste0(file, ".110-500"))
plotit(log[500:1000,], outpref=paste0(file, ".500-1000"))
