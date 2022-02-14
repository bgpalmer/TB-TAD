

###############
###############
##
## Extract the data for these samples
##
###############
###############

bdsGR <-bdsGR_GM12878 #select the cell line you want

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