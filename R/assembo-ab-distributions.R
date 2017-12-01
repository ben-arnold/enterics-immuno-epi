#-----------------------------
# assembo-ab-distributions.R
#
# Plot antibody distributions
# from the Assembo cohort
# using histograms and 
# smoothed kernel densities
#-----------------------------


#-----------------------------
# preamble
#-----------------------------
library(tidyverse)


#-----------------------------
# load the data
#-----------------------------

d <- readRDS("~/dropbox/kisumu-sero/data/final/assembo_sero.rds")

mbavars <- c("vsp3","vsp5","cp17","cp23","leca","rotavirus","p18","p39","etec","cholera")
mbalabs <- c("Giardia VSP-3","Giardia VSP-5","Cryptosporidium Cp17","Cryptosporidium Cp23",
             "E. histolytica LecA","Rotavirus","Campylobacter p18","Campylobacter p39","ETEC toxin beta subunit","Cholera toxin beta subunit")
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
# add cutoff values for
# those with an ROC cutoff
#-----------------------------
dl$roccut <- NA
dl$roccut[dl$antigen %in% "vsp3"] <- 542
dl$roccut[dl$antigen %in% "vsp5"] <- 421
dl$roccut[dl$antigen %in% "cp17"] <- 375
dl$roccut[dl$antigen %in% "cp23"] <- 1336

#-----------------------------
# group antigens by pathogen
# for plot aesthetics
# label antigens with fancier names
#-----------------------------
dl$pathogen <- NA
dl$pathogen[dl$antigen %in% c("vsp3","vsp5")] <- "Giardia"
dl$pathogen[dl$antigen %in% c("cp17","cp23")] <- "Cryptosporidium"
dl$pathogen[dl$antigen %in% c("leca")] <- "E. histolytica"
dl$pathogen[dl$antigen %in% c("rotavirus")] <- "Rotavirus"
dl$pathogen[dl$antigen %in% c("p18","p39")] <- "Campylobacter"
dl$pathogen[dl$antigen %in% c("etec")] <- "ETEC"
dl$pathogen[dl$antigen %in% c("cholera")] <- "V. cholerae"
dl <- dl %>%
  mutate(pathogen = factor(pathogen,levels=c("Giardia","Cryptosporidium","E. histolytica","Rotavirus","Campylobacter","ETEC","V. cholerae")),
         antigen=factor(antigen,levels=mbavars,labels=mbalabs))

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
p <- ggplot(data=dl,aes(x=log10mfi,group=pathogen,color=pathogen,fill=pathogen)) +
  facet_wrap(~antigen,nrow=5,ncol=2,scales="free_y") +
  geom_histogram(aes(y=..density..),bins=50,alpha=0.7) +
  geom_density(aes(fill=NULL),color="black") +
  geom_vline(aes(xintercept=log10(roccut))) +
  labs(x="Log10 Luminex Response (MFI-bg)") +
  scale_fill_manual(values=pcols) +
  scale_color_manual(values=pcols) +
  theme_minimal() +
  theme(legend.position = "none")

p

ggsave("~/dropbox/assembo/results/figs/assembo-ab-distributions.pdf",plot=p,device="pdf",width=7,height=14)



