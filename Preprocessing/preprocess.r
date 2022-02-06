# this file is for preprocessing the data for

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)

#######################################
#####   methylation infos
#######################################
# library(FDb.InfiniumMethylation.hg19)

prob_data<-read.table(paste0(Sys.getenv('DATADIR'), '/Methylation/GSE62111_series_matrix.txt'),header=T,sep = '\t',stringsAsFactors = F, fill=TRUE, comment.char='!')

# library(FDb.InfiniumMethylation.hg19)
# hm450 <- get450k()
hm450.filepath <- paste0(Sys.getenv('DATADIR'), '/Methylation/hm450.Rdata')
# save(hm450,file=hm450.filepath)
load(hm450.filepath)

head(hm450)

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


#######################################
#####   Save the file
#######################################

# library(tidyverse)

# tb <- tibble(
# 	date = as.Date(date_time),
# 	hour = hour(date_time),
# 	minute = minute(date_time),
# 	second = second(date_time)
# ) %>% 
# mutate(
# 	format_date = format(date, "%m/%d/%Y"),
# 	format_hour = paste(hour, minute, second, sep = ":")
# )

save(list=ls(),file='pre_info_08132020.Rdata')