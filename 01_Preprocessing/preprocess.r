# this file is for preprocessing the data for

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)

#Get TAD boundaries from TAD file
# TAD<- read_table(paste0(Sys.getenv('DATADIR'), "/Hi-C/GSM455133_30E0LAAXX.1.maq.hic.summary.binned.txt"),col_names = F)
TAD<- read_table(paste0(Sys.getenv('DATADIR'), "/TADs/GM12878_Lieberman-raw_TADs.txt"), col_names = F)
# head(TAD)
# nrow(TAD)
#colnames(TAD)[1:3]<- c("X1","X2","X3") #chr, start, stop
bdsdf <-matrix("0",nrow=nrow(TAD),ncol=3)
colnames(bdsdf)<-c("seqnames","start","end")
for(i in 1:(nrow(TAD)-1)){
  if(TAD[i,"X1"] == TAD[i+1,"X1"])
    {
    bdsdf[i,"seqnames"]<-as.character(TAD[i,"X1"])
    bdsdf[i,"start"]<-(as.numeric(TAD[i,"X3"])+as.numeric(TAD[i+1,"X2"])+1)/2 - 100000 #boundary width is 200k
    bdsdf[i,"end"]<-(as.numeric(TAD[i,"X3"])+as.numeric(TAD[i+1,"X2"])+1)/2 + 100000
    }
}

bdsdf<-bdsdf[-(which(bdsdf[,"seqnames"]=="0")),] #remove the rows inbetween chromosomes
bdsGR <- GRanges(seqnames = bdsdf[,"seqnames"],
                      ranges = IRanges(start=as.numeric(bdsdf[,"start"]),end=as.numeric(bdsdf[,"end"])),
                      strand='*')
                      
bdsGR_GM12878 <- bdsGR

# save(list=ls(),file='bdsdf.Rdata')



#######################################
#####   methylation infos
#######################################
library(FDb.InfiniumMethylation.hg19)

prob_data<-read.table(paste0(Sys.getenv('DATADIR'), '/Methylation/GSE62111_series_matrix.txt'),header=T,sep = '\t',stringsAsFactors = F, fill=TRUE, comment.char='!')

library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
hm450.filepath <- paste0(Sys.getenv('DATADIR'), '/Methylation/hm450.Rdata')
save(hm450,file=hm450.filepath)
load(hm450.filepath)

# head(hm450)

y<-as.data.table(hm450)
y[,'ID_REF']<-names(hm450)

colnames(prob_data)
tmp<-prob_data[,c("ID_REF","GSM1519791")] #GM12878 methylation https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62111

names(tmp)<-c('ID_REF','v1')
prob_data_all<-left_join(tmp,y=y)

probGR<-GRanges(seqnames=(prob_data_all$seqnames)
                ,ranges=IRanges(start=prob_data_all$start,end=prob_data_all$start)
                ,strand='*'
                ,v1=prob_data_all$v1
                ,p1=0)
me <- probGR

chrnames<-c(paste('chr',1:22,sep=''),'chrX')



##get center and chr_size information for each chrs
chr_info_1<-read_tsv(paste0(Sys.getenv('DATADIR'),'/Sizes/hg19_chrom_centros.txt'), col_names = F)
chr_info_2<-read_tsv(paste0(Sys.getenv('DATADIR'),'/Sizes/hg19_chrom_sizes.txt'),col_names = F)


# #######################################
# #####   Save the file
# #######################################

# # library(tidyverse)

# # tb <- tibble(
# # 	date = as.Date(date_time),
# # 	hour = hour(date_time),
# # 	minute = minute(date_time),
# # 	second = second(date_time)
# # ) %>% 
# # mutate(
# # 	format_date = format(date, "%m/%d/%Y"),
# # 	format_hour = paste(hour, minute, second, sep = ":")
# # )

save(list=ls(),file='pre_info_08132020.Rdata')