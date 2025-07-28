#==get discards
ret_cat<-read.csv("data/new_catch/retained_catch.csv")
ret_cat$fish<-substring(ret_cat$fishery,1,2)
tot_cat<-read.csv("data/new_catch/total_catch.csv")
disc<-read.csv("data/new_catch/bssc_discards.csv")

library(dplyr)
library(ggplot2)
library(reshape2)
library(png)
library(grid)
library(tidyr)

fems<-filter(tot_cat,group=='female')
fems$fish<-substring(fems$fishery,1,2)
  
ggplot(fems)+
  geom_line(aes(x=crab_year,y=total_catch_wt,col=fish,group=fish))+
  theme_bw()

males<-filter(tot_cat,group!='female')%>%
  group_by(crab_year,fishery)%>%
  summarize(tot_cat=sum(total_catch_wt)/1000)
males$fish<-substring(males$fishery,1,2)

png("plots/male_discards.png")
ggplot(males)+
  geom_line(aes(x=crab_year,y=tot_cat,col=fish,group=fish))+
  theme_bw()+ylab("Total catch")
dev.off()

ggplot(filter(males,fish!="QO"))+
  geom_line(aes(x=crab_year,y=tot_cat,col=fish,group=fish))+
  theme_bw()+ylab("Total catch")

directed_discard<-merge(ret_cat,filter(males,fish=="QO"))
directed_discard$tot_retained_wt <-directed_discard$tot_retained_wt /1000
directed_discard$subtract_disc<-directed_discard$tot_cat -directed_discard$tot_retained_wt  


chrp<-filter(directed_discard,subtract_disc>0)
hist(chrp$subtract_disc/chrp$tot_cat)
med_disc_rate<-median(chrp$subtract_disc/chrp$tot_cat)

directed_discard$use_disc<-directed_discard$subtract_disc
directed_discard$use_disc[which(directed_discard$subtract_disc<0)]<-directed_discard$tot_retained_wt[which(directed_discard$subtract_disc<0)]*med_disc_rate


write.csv(cbind(directed_discard$crab_year,directed_discard$use_disc,filter(fems,fish=='QO')$total_catch_wt/1000),
          "data/derived/dir_disc_m_f.csv")

#=============================
# directed size comps
#==============================
ret_cat_sc<-read.csv("data/new_catch/retained_catch_composition.csv")
tot_cat_sc<-read.csv("data/new_catch/directed_total_composition.csv")

#==retained
data<-ret_cat_sc%>%
  group_by(crab_year,size)%>%
  summarize(tot_crab=sum(total))

# Define midpoints
midpoints <- seq(27.5, 132.5, by = 5)

# Calculate bin edges from midpoints
bin_width <- 5
bin_edges <- c(midpoints - bin_width / 2, last(midpoints) + bin_width / 2)

# Create bins
data <- data %>%
  mutate(bin = cut(size, breaks = bin_edges, include.lowest = TRUE, right = TRUE,labels=midpoints))

# Summarize within bins
binned_data <- data %>%
  group_by(crab_year, bin) %>%
  summarize(tot_crab_sum = sum(tot_crab), .groups = 'drop')

# Normalize the data within each year
normalized_data_ret <- binned_data %>%
  group_by(crab_year) %>%
  mutate(total_crab_year = sum(tot_crab_sum),
         normalized_crab = tot_crab_sum / total_crab_year) %>%
  ungroup()


output<-dcast(normalized_data_ret,crab_year~bin,id.var='normalized_crab')
output[is.na(output)]<-0
rownames(output)<-output[,1]
write.csv(output[,-1],"data/derived/ret_sc.csv")

#==retained
data<-tot_cat_sc%>%
  group_by(crab_year,size)%>%
  summarize(tot_crab=sum(total))

# Create bins
data <- data %>%
  mutate(bin = cut(size, breaks = bin_edges, include.lowest = TRUE, right = TRUE,labels=midpoints))

# Summarize within bins
binned_data <- data %>%
  group_by(crab_year, bin) %>%
  summarize(tot_crab_sum = sum(tot_crab), .groups = 'drop')

# Normalize the data within each year
normalized_data_tot <- binned_data %>%
  group_by(crab_year) %>%
  mutate(total_crab_year = sum(tot_crab_sum),
         normalized_crab = tot_crab_sum / total_crab_year) %>%
  ungroup()

output<-dcast(normalized_data_tot,crab_year~bin,id.var='normalized_crab')
output[is.na(output)]<-0
rownames(output)<-output[,1]
write.csv(output[,-1],"data/derived/tot_sc.csv")













#=======================================================================
# bycatch length composition
# these data come from the "Observer data" tab on AKfin
# Enter 'NORPAC Length Report - Haul & Length"
# Download the data for snow crab
# Time period should be July 1 (last year) to June 30 (this year)
#=======================================================================

LenDatBig<-read.csv("data/norpac_length_report/norpac_length_report.csv",skip=6)
LenDatBig$Haul.Offload.Date<-strptime(LenDatBig$Haul.Offload.Date,format="%d-%b-%y")
range(LenDatBig$Haul.Offload.Date)

#==CHECK DATES, EH!
LenDat<-LenDatBig[LenDatBig$Haul.Offload.Date >= "2023-07-01" & LenDatBig$Haul.Offload.Date <= "2024-06-30" & LenDatBig$Species.Name=="OPILIO TANNER CRAB",]

#==males old shell discard
LengthBins			<-seq(25,135,5)
BycatchFem			<-rep(0,length(LengthBins))
BycatchMale			<-rep(0,length(LengthBins))

for(y in 1:(length(LengthBins)-1))
{
  BycatchFem[y]	<-sum(LenDat$Frequency[LenDat$Length..cm.>=LengthBins[y] & LenDat$Length..cm.<LengthBins[y+1] & LenDat$Sex=="F"])
  BycatchMale[y]	<-sum(LenDat$Frequency[LenDat$Length..cm.>=LengthBins[y] & LenDat$Length..cm.<LengthBins[y+1] & LenDat$Sex=="M"])
}

upperBnd							<-130
BycatchFem[which(LengthBins==upperBnd)]	<-sum(BycatchFem[(which(LengthBins==upperBnd)):length(LengthBins)])
BycatchFem						<-BycatchFem[1:which(LengthBins==upperBnd)]
BycatchMale[which(LengthBins==upperBnd)]	<-sum(BycatchMale[(which(LengthBins==upperBnd)):length(LengthBins)])
BycatchMale						<-BycatchMale[1:which(LengthBins==upperBnd)]

par(mfrow=c(1,2))
barplot(BycatchFem,names.arg=LengthBins,xlab="Carapace width (mm)",ylab="Count")
barplot(BycatchMale)

write.table(rbind(BycatchFem,BycatchMale),"data/derived/bycatch_len_comps_f_then_m.txt",
            row.names=FALSE,col.names=F)

#=======================================================================
# bycatch numbers
# these data come from the "Observer data" tab on AKfin
# Enter 'NORPAC catch Report"
# Download the data for snow crab in BS of BSAI
# Time period should be July 1 (last year) to June 30 (this year)
#=======================================================================


#==if this doesn't work, check the number of lines to skip
bycatch_dat<-read.csv("data/norpac_catch_report/norpac_catch_report.csv",skip=7)
temp<-strptime(bycatch_dat$Haul.Date,format="%d-%b-%y")
bycatch_dat$Haul.Date<-substr(temp,start=1,stop=10)

#==bycatch numbers total
bycatchDat<-bycatch_dat[bycatch_dat$Haul.Date >= "2023-07-01" & bycatch_dat$Haul.Date <= "2024-06-30"& 
                          bycatch_dat$Species.Name=="OPILIO TANNER CRAB" & bycatch_dat$Gear.Description!="POT OR TRAP" ,]
bycatch_num_tot<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
bycatch_wt_tot<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)


#=========================================
# all years bycatch
#==================================

bycatch_dat_big<-bycatch_dat
#temp<-strptime(bycatch_dat_big$Haul.Date,format="%y-%b-%d")
#bycatch_dat_big$Haul.Date<-substr(temp,start=1,stop=10)
bycatch_dat_big$crab.year<-bycatch_dat_big$Year

#==bycatch numbers total
bycatch_year    <-sort(unique(bycatch_dat_big$Year))
bycatch_num_tot		<-rep(0,length(bycatch_year))
bycatch_wt_tot		<-rep(0,length(bycatch_year))

for(y in 1:(length(bycatch_year)-1))
{  
  bycatchDat<-bycatch_dat_big[bycatch_dat_big$Haul.Date >= paste(bycatch_year[y],"-07-01",sep="") & bycatch_dat_big$Haul.Date <= paste(bycatch_year[y]+1,"-06-30",sep="") & 
                                bycatch_dat_big$Species.Name=="OPILIO TANNER CRAB" ,]
  bycatch_num_tot[y]<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  bycatch_wt_tot[y]<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)  
}



#===add in all of the crab fisheries
all_crab_bycatch<-filter(males,fish!="QO")%>%
  group_by(crab_year)%>%
  summarize(all_crab_bycatch=sum(tot_cat,na.rm=T)*0.3)

#==need to apply .3 mortality to pot and .8 to trawl
nondir<-cbind(bycatch_year,bycatch_wt_tot)
nondir[,2]<-(nondir[,2]/1000000)*0.8
colnames(nondir)[1]<-'crab_year'

all_nondir_bycatch<-merge(all_crab_bycatch,nondir)
all_nondir_bycatch$tot_nondir<-all_nondir_bycatch$all_crab_bycatch+all_nondir_bycatch$bycatch_wt_tot
write.csv(all_nondir_bycatch,"data/derived/bycatch_wt_total.csv")


#==bycatch by gear type
library(dplyr)
library(ggplot2)
for(y in 1:(length(bycatch_year)-1))
  bycatch_dat_big$crab.year[bycatch_dat_big$Haul.Date >= paste(bycatch_year[y],"-07-01",sep="") & bycatch_dat_big$Haul.Date <= paste(bycatch_year[y]+1,"-06-30",sep="")]<-bycatch_year[y]

in_dat<-bycatch_dat_big[,-24]
temp<- in_dat %>%
  group_by(crab.year,Gear.Description) %>%
  summarise(Bycatch = sum(Extrapolated.Number))

temp$Year<-as.numeric(temp$Year)

png("plots/bycatch.png",height=8,width=8,res=400,units='in')
ggplot(temp)+
  geom_line(aes(x=crab.year,y=Bycatch,col=Gear.Description))+
  geom_point(aes(x=crab.year,y=Bycatch,col=Gear.Description))+
  theme_bw()+theme(legend.position=c(.8,.8))
dev.off()



#================================
# Growth data

grow<-read.csv("data/SnowCrabGrowthMaster.csv")
use_grow<-filter(grow,Legs_missing_premolt==0)
use_grow$molt_inc<-use_grow$Postmolt_CW-use_grow$Premolt_CW
new_grow<-cbind(use_grow[,c(6,3,17)],rep(0.03,nrow(use_grow)))
write.csv(new_grow,'data/newgrowth.csv')
ass_grow<-read.csv('data/growth_ass.csv')
keepers<-ass_grow
colnames(keepers)<-c("pre","sex","inc","cv")
colnames(new_grow)<-c("pre","sex","inc","cv")
for(x in 1:nrow(new_grow))
{
  check<-new_grow[x,]
  takeit<-1
  for(y in 1:length(ass_grow))
  {
    matched<-sum(!is.na(match(check,ass_grow[y,])))
    if(matched==4)
      takeit<-0
  }  
  if(takeit==1)
    keepers<-rbind(keepers,check)
}

write.csv(keepers,'data/newgrowth2.csv')

