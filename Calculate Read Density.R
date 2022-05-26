# Read diversity
rm(list=ls())
library(dplyr)
Root <- read.delim("~/Desktop/BACKUP/Root_P_sites_all_4Col.bed",header=F,sep="\t",stringsAsFactors = F,skip=0) #2,157,859 rows, sum total reads 132,308,885
Shoot <- read.delim("~/Desktop/BACKUP/Shoot_P_sites_all_4Col.bed",header=F,sep="\t",stringsAsFactors = F,skip=0) #1,908,790 rows, sum total reads 106,752,951
NewCTRL <- read.delim("/Users/wu/Desktop/CTRL_v1/CTRL_TPM0.25_P_sites_sort_count",header=F,sep="\t",stringsAsFactors = F,skip=0) #11,209,299 rows, sum total reads 298,004,424
NEB <- read.delim("~/Desktop/New_Riboseq/NEB123/P_sites_all_processed",header=F,sep="\t",stringsAsFactors = F,skip=0) # rows, sum total reads 


RootChr1 <- Root %>% filter(V2 %in% c(1,2,3))
nrow(RootChr1) #1298078
ShootChr1 <- Shoot %>% filter(V2 %in% c(1,2,3))
nrow(ShootChr1) #1143017
NewCTRLChr1 <- NewCTRL %>% filter(V2 %in% c(1,2,3))
nrow(NewCTRLChr1) #6772778
NEBChr1 <- NEB %>% filter(V2 %in% c(1,2,3))
nrow(NEBChr1) #3272819

RootChr1R <- as.data.frame(lapply(RootChr1, rep, RootChr1$V1))
nrow(RootChr1R) #75877138
ShootChr1R <- as.data.frame(lapply(ShootChr1, rep, ShootChr1$V1))
nrow(ShootChr1R) #61771627
NewCTRLChr1R <- as.data.frame(lapply(NewCTRLChr1, rep, NewCTRLChr1$V1))
nrow(NewCTRLChr1R) #182677883
NEBChr1R <- as.data.frame(lapply(NEBChr1, rep, NEBChr1$V1))
nrow(NEBChr1R) #21693744

NUM <- seq(1e+6,2e+7, by=1e+6)

RootChr1Cdf <- c()
for(i in 1:20) {
  print(i)
  RootChr1C <- RootChr1R[sample(nrow(RootChr1R), NUM[i]), ]
  res <- RootChr1C[!duplicated(RootChr1C), ]
  RootChr1Cdf[i] <- nrow(res)
}
RootChr1Cdf

ShootChr1Cdf <- c()
for(i in 1:20) {
  print(i)
  ShootChr1C <- ShootChr1R[sample(nrow(ShootChr1R), NUM[i]), ]
  res <- ShootChr1C[!duplicated(ShootChr1C), ]
  ShootChr1Cdf[i] <- nrow(res)
}
ShootChr1Cdf

NewCTRLChr1Cdf <- c()
for(i in 1:20) {
  print(i)
  NewCTRLChr1C <- NewCTRLChr1R[sample(nrow(NewCTRLChr1R), NUM[i]), ]
  res <- NewCTRLChr1C[!duplicated(NewCTRLChr1C), ]
  NewCTRLChr1Cdf[i] <- nrow(res)
}
NewCTRLChr1Cdf

NEBChr1Cdf <- c()
for(i in 1:20) {
  print(i)
  NEBChr1C <- NEBChr1R[sample(nrow(NEBChr1R), NUM[i]), ]
  res <- NEBChr1C[!duplicated(NEBChr1C), ]
  NEBChr1Cdf[i] <- nrow(res)
}
NEBChr1Cdf

dfNEB <- data.frame(total_read=NUM, distinct_p_site=NEBChr1Cdf,Samples="At_NEB")
dfCTRL <- data.frame(total_read=NUM, distinct_p_site=NewCTRLChr1Cdf,Samples="At_CTRL")
dfRoot <- data.frame(total_read=NUM, distinct_p_site=RootChr1Cdf,Samples="At_Root")
dfShoot <- data.frame(total_read=NUM, distinct_p_site=ShootChr1Cdf,Samples="At_Shoot")

df <- rbind(dfNEB,dfCTRL,dfRoot,dfShoot)
save(df,"~/Desktop/readDiversity.RData")
df$distinct_p_site <- df$distinct_p_site/1000000
p3 <- ggplot(df, aes(x=total_read,y=distinct_p_site, colour = Samples)) +
  geom_line(aes(color=Samples),size=0.7) +
  geom_point(aes(color=Samples)) +
  theme_bw() +
  labs(x="Million footprint counts",y="Million distinct P-sites in 20M samples")
p3
ggsave("~/Desktop/Footprint_diversity_NEB.pdf") 
