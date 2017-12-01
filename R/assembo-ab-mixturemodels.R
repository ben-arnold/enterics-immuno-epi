#-----------------------------
# assembo-ab-mixturemodels.R
#
# Determine seropositivity
# using gaussian mixture models
#-----------------------------


#-----------------------------
# preamble
#-----------------------------
rm(list=ls())
library(tidyverse)
library(mixtools)

# set up for parallel computing
# configure for a laptop (use only 3 cores)
library(foreach)
library(doParallel)
registerDoParallel(cores=3)


#-----------------------------
# load the data
#-----------------------------

d <- readRDS("~/dropbox/assembo/data/final/assembo_sero.rds")

mbavars <- c("vsp3","vsp5","cp17","cp23","leca","rotavirus","p18","p39","etec","cholera")
mbalabs <- c("Giardia VSP-3","Giardia VSP-5","Cryptosporidium Cp17","Cryptosporidium Cp23",
             "E. histolytica LecA","Rotavirus","Campylobacter p18","Campylobacter p39","ETEC toxin beta subunit","Cholera toxin beta subunit")

#-----------------------------
# restrict to measurements <= 12 months old
#-----------------------------
# d <- d %>%
#   filter(age<=12)

#-----------------------------
# reshape long and convert
# to log10 values
#-----------------------------
dl <- d %>%
  select(childid,time,sex,age,tr,mbavars) %>%
  gather(antigen,mfi,-childid,-time,-sex,-age,-tr)

# set negative and zero values to 1 before the log10 transform
dl <- dl %>%
  mutate(mfi = ifelse(mfi<=0,1,mfi),
         log10mfi = log10(mfi),
         antigen=factor(antigen,levels=mbavars)
  )
         
#-----------------------------
# mixture models
#-----------------------------
fitmix <- function(x,lambda,k) {
  mixfit <- normalmixEM(x,lambda=lambda,k=k)
  mixcut <- mixfit$mu[order(mixfit$mu)][1]+3*mixfit$sigma[order(mixfit$mu)][1]
  cat(summary(mixfit),"\nCutoff value:",mixcut,"(log10 MFI), or:", round(10^mixcut),"(MFI)\n\n")
  # pull out fitted densities
  denmat <- matrix(NA,nrow=length(x),ncol=k)
  for(i in 1:k) {
    denmat[,i] <- mixfit$lambda[i] * dnorm(x,mean=mixfit$mu[i],sd=mixfit$sigma[i])
  }
  denmat <- data.frame(denmat)
  colnames(denmat) <- paste("den",1:k,sep="")
  # return original values plus fitted densities in a dataframe 
  # also return the cutoff value and normalmixEM object
  xden <- data.frame(x=x,denmat)
  list(xden=xden,mixcut=mixcut,mixfit=mixfit)
  
}

## mixture model fits -- store densities of mixture distributions
mixdens <- foreach(ab=mbavars,.combine=rbind) %dopar% {
  mixf <- fitmix(x=dl$log10mfi[dl$antigen==ab],lambda=0.5,k=2)
  di <- mixf$xden
  di$antigen=ab
  di <- di %>% arrange(x)
  di
}

## mixture model fits -- store cutoff values
mixcuts <- foreach(ab=mbavars,.combine=rbind) %do% {
  mixf <- fitmix(x=dl$log10mfi[dl$antigen==ab],lambda=0.5,k=2)
  di <- data.frame(antigen=ab,mixcut=mixf$mixcut)
  di
}

# Distributions for p18,ETEC and Cholera do not allow for reasonable mixture cutoffs
mixcuts$mixcut[mixcuts$antigen %in% c("p18","etec","cholera")] <- NA

#-----------------------------
# merge in the mixture cutoffs
#-----------------------------
dl <- left_join(dl,mixcuts,by=c("antigen"))

#-----------------------------
# add cutoff values for
# those with an ROC cutoff
#-----------------------------
dl <- dl %>%
  mutate(roccut = NA)
dl$roccut[dl$antigen %in% "vsp3"] <- log10(542)
dl$roccut[dl$antigen %in% "vsp5"] <- log10(421)
dl$roccut[dl$antigen %in% "cp17"] <- log10(375)
dl$roccut[dl$antigen %in% "cp23"] <- log10(1336)

#-----------------------------
# save the results for later
# use in calculating seroprevalence
# and seroincidence
#-----------------------------
save.image(file="~/dropbox/assembo/results/raw/assembo-ab-mixturemodels.RData")

#-----------------------------
# group antigens by pathogen
# for plot aesthetics
# label antigens with fancier names
#-----------------------------

dplot <- dl

dplot$pathogen <- NA
dplot$pathogen[dplot$antigen %in% c("vsp3","vsp5")] <- "Giardia"
dplot$pathogen[dplot$antigen %in% c("cp17","cp23")] <- "Cryptosporidium"
dplot$pathogen[dplot$antigen %in% c("leca")] <- "E. histolytica"
dplot$pathogen[dplot$antigen %in% c("rotavirus")] <- "Rotavirus"
dplot$pathogen[dplot$antigen %in% c("p18","p39")] <- "Campylobacter"
dplot$pathogen[dplot$antigen %in% c("etec")] <- "ETEC"
dplot$pathogen[dplot$antigen %in% c("cholera")] <- "V. cholerae"
dplot <- dplot %>%
  mutate(pathogen = factor(pathogen,levels=c("Giardia","Cryptosporidium","E. histolytica","Rotavirus","Campylobacter","ETEC","V. cholerae")),
         antigen=factor(antigen,levels=mbavars,labels=mbalabs))


# for the mixture distributions, set p18, ETEC, cholera to missing
plotmixdens <- mixdens %>%
  mutate(den1p = ifelse(antigen %in% c("p18","etec","cholera"),NA,den1),
         den2p = ifelse(antigen %in% c("p18","etec","cholera"),NA,den2)
         ) %>%
  mutate(antigen=factor(antigen,levels=mbavars,labels=mbalabs),
         pathogen="Giardia"
         )

#-----------------------------
# plot distributions
#-----------------------------

# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"

pcols <- c(cred,corange,cchartr,cgreen,cteal,cblue,cmagent)

# group=antigen,color=antigen,fill=antigen
p <- ggplot(data=dplot,aes(x=log10mfi,group=pathogen,color=pathogen,fill=pathogen)) +
  facet_wrap(~antigen,nrow=5,ncol=2,scales="free_y") +
  # plot empirical distribution and smooth
  geom_histogram(aes(y=..density..),bins=50,alpha=0.7) +
  geom_density(aes(fill=NULL),color="black") +
  # plot fitted mixture distributions
  geom_line(data=plotmixdens,aes(x=x,y=den1p),color=cgrey) +
  geom_line(data=plotmixdens,aes(x=x,y=den2p),color=cgrey) +
  # add vertical lines for the cutoff
  geom_vline(aes(xintercept=roccut)) +
  geom_vline(aes(xintercept=mixcut),linetype=2) +
  # labels and formatting
  scale_x_continuous(limits = c(0,4.5),breaks = 0:4,labels=0:4) +
  coord_cartesian(ylim=c(0,2)) +
  labs(x="Log10 Luminex Response (MFI-bg)") +
  scale_fill_manual(values=pcols) +
  scale_color_manual(values=pcols) +
  theme_minimal() +
  theme(legend.position = "none")

p

ggsave("~/dropbox/assembo/results/figs/assembo-ab-mixture-distributions.pdf",plot=p,device="pdf",width=7,height=14)











