

###############
###############
##
## Extract the data for these samples
##
###############
###############

library(GenomicRanges)

load('../01_Preprocessing/pre_info_08132020.Rdata') #preprocessed ChIPseq and epigenomic and genomic data


region_bin_len=10000 #width of each genomic region
test_ratio=3/10 #ratio of the testing set
neighbors=10 #number of neghboring bins to consider

region_GR<-GRanges()
for (kk in 1:23){
  chrklen<-chr_info_2$X2[chr_info_2$X1%in% chrnames[kk]]
  chrkseq<-seq(from=1,by=region_bin_len,to=chrklen)
  bin_num<-length(chrkseq)
  tr<-rep(1,(bin_num-1))
  tr[sample(1:(bin_num-1),size = test_ratio*(bin_num-1))]<-0
  tmpGR<-GRanges(seqnames = chrnames[kk]
                 ,ranges = IRanges(start=chrkseq[1:(bin_num-1)],end = chrkseq[2:bin_num])
                 ,strand = '*'
                 ,train= tr
  )
  region_GR<-append(region_GR,tmpGR)
}

chr_info_tc <- chr_info_1[which(chr_info_1$type %in% c("centromere","telomere")),] # identify telomere and centromere regions

print(chr_info_tc)

chr_info_tc <- GRanges(seqnames = chr_info_tc$chrom,
                       ranges = IRanges(start=chr_info_tc$chromStart, end = chr_info_tc$chromEnd),
                       strand = "*",
                       type = chr_info_tc$type)

tmp_ov <- findOverlaps(region_GR, chr_info_tc,type="any")
region_GR$telocentro <- "0"
region_GR$telocentro[tmp_ov@from]<- "1"


bdsGR <-pre_info_08132020.Rdata #select the cell line you want

gen_data_1<-function(train_region=region_GR,by_len=10000,n_neighbors=10){
  ana_list<-list()
  for(kk in 1:23){
    
    #construct data for modeling
    chkk<-chrnames[kk]
    bdsGR_chrkk<-bdsGR[which(bdsGR@seqnames == chrnames[kk])]
    ana_data<-data.frame(start=seq(from=0,to=max(bdsGR@ranges@start+bdsGR@ranges@width),by=by_len)+1
                         ,end=seq(from=by_len,to=max(bdsGR@ranges@start+bdsGR@ranges@width)+by_len,by=by_len))
    
    #using GRanges ****************************************************************
    anGR<-GRanges(seqnames=chkk
                  ,ranges=IRanges(start=ana_data$start,end=ana_data$end)
                  ,strand='*')
    
    tmp<-findOverlaps(anGR,train_region,type='within',ignore.strand=T)
    anGR<-anGR[tmp@from]


	#################################################
	###     add methylation values
	#################################################
	tmp<-findOverlaps(me,anGR,type='within',ignore.strand=T)
	prob_valid <- me$p1<0.05 #p-value
	prob_v<-cbind(me$v1) #methylation beta value
	tmp1<-data.table(prob=(rowSums(prob_v*prob_valid)/prob_valid)[tmp@from],id=tmp@to)
	mcols(anGR)[,'me']<-NA
	tmp2<-tmp1[,.(mv=mean(prob,na.rm=T)),by=id]
	mcols(anGR)[,'me'][tmp2$id]<-tmp2$mv


all_data_1 <- gen_data_1(region_GR,by_len=region_bin_len,n_neighbors=neighbors)
save(all_data_1,region_GR,file="all_data_1_GM12878.Rdata")