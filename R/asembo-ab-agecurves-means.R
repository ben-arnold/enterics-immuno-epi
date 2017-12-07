#-----------------------------
# asembo-ab-agecurves-means.R
#
# Estimate age-dependent antibody
# curves for log10 antibody levels
#-----------------------------

#-----------------------------
# preamble
#-----------------------------
library(tidyverse)

# set up for parallel computing
# configure for a laptop (use only 3 cores)
library(foreach)
library(doParallel)
registerDoParallel(cores=3)



#-----------------------------
# load the data
#-----------------------------

d <- readRDS("~/dropbox/asembo/data/final/asembo_sero.rds")

d <- d %>%
  ungroup() %>%
  mutate(childid=factor(childid))

# filter out 6 measurements with no age for now
d <- d %>%
  filter(!is.na(age))


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
         antigenf=factor(antigen,levels=mbavars,labels=mbalabs))

#-----------------------------
# estimate curves with 
# splines
#-----------------------------


#----------------------------------
# simulataneous CIs for GAMs
# estimated by resampling the 
# Baysian posterior estimates of
# the variance-covariance matrix
# assuming that it is multivariate normal
# much better coverage than pointwise
# see: http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#----------------------------------

gamCI <- function(m,newdata,nreps=10000) {
  require(mgcv)
  require(dplyr)
  Vb <- vcov(m,unconditional = TRUE)
  pred <- predict(m, newdata, se.fit = TRUE)
  fit <- pred$fit
  se.fit <- pred$se.fit
  BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg <- predict(m, newdata, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- data.frame(newdata,fit=pred$fit,se.fit=pred$se.fit)
  pred <- mutate(pred,
                 uprP = fit + (2 * se.fit),
                 lwrP = fit - (2 * se.fit),
                 uprS = fit + (crit * se.fit),
                 lwrS = fit - (crit * se.fit)
  )
  return(pred)
}


# fit cubic spline over age, separately for each antigen
library(mgcv)
gamfits <- foreach(ab=levels(dl$antigen),.combine=rbind) %dopar% {
    pd <- filter(dl, antigen==ab)
    gfit <- gam(log10mfi~s(age, bs="cr"),data=pd)
    gsci <- gamCI(m=gfit,newdata=pd,nreps=1000)
    gsci$antigen <- ab
    return(gsci)
  }
  

#-----------------------------
# estimate curves with 
# ensemble machine learning
# including a library with
# GLM, GAM, Loess, and 
# extreme gradient boosting
#-----------------------------
detach("package:mgcv", unload=TRUE)
library(tmleAb)
set.seed(12345)
slfits <-  foreach(ab=levels(dl$antigen),.combine=rbind) %dopar% {
  pd <- filter(dl, antigen==ab)
  slfit <- agecurveAb(Y=pd$log10mfi,
                      Age=pd$age,
                      id=pd$childid,
                      SL.library = c("SL.glm","SL.gam","SL.loess","SL.xgboost"))
  sld <- data.frame(age=slfit$Age,fit=slfit$pY)
  sld <- sld %>% group_by(age) %>% filter(row_number()==1) %>% ungroup()
  slres <- left_join(pd,sld,by=c("age"))
  return(slres)
}


#-----------------------------
# plot results
#-----------------------------

log10labs <- c( 
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)
)

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

p <- ggplot(data=gamfits,aes(x=age,group=pathogen,color=pathogen,fill=pathogen)) +
  # approximate simultaneous CI from spline fit
  geom_ribbon(aes(ymin=lwrS,ymax=uprS),alpha=0.2,color=NA) +
  # actual measurements
  geom_jitter(aes(y=log10mfi),width=0.2,height=0,color="black",alpha=0.4,size=0.3) +
  # ensemble fit
  geom_line(data=slfits,aes(x=age,y=fit),color=cgrey) +
  # spline fit
  geom_line(aes(y=fit)) +
  # plot aesthetics
  facet_wrap(~antigenf,nrow=5,ncol=2) +
  scale_x_continuous(breaks=seq(4,18,by=2))+
  scale_y_continuous(breaks=0:4,labels=log10labs)+
  coord_cartesian(ylim=c(0,4.5))+
  labs(x="Age in months",y="Luminex Response (MFI-bg)") +
  scale_fill_manual(values=pcols) +
  scale_color_manual(values=pcols) +
  theme_minimal() +
  theme(legend.position = "none")

p

ggsave("~/dropbox/asembo/results/figs/asembo-ab-agecurves-mean.pdf",plot=p,device="pdf",width=7,height=14)



