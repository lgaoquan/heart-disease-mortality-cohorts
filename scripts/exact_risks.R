suppressPackageStartupMessages({library(survival);library(data.table)})
root<-Sys.getenv('ANALYSIS_OUTPUT_ROOT', unset='analysis_output');b<-fread(file.path(root,'private/baseline.csv'));b[,duration:=round(duration,8)];curves<-list()
for(co in c('CHARLS','ELSA','HRS','SHARE')){
 x<-as.data.frame(b[cohort==co & complete==TRUE]);x$baseline_heart<-factor(x$baseline_heart,levels=c(0,1))
 sf<-survfit(Surv(duration,event)~baseline_heart,data=x,conf.type='log-log');sm<-summary(sf,censored=TRUE)
 curves[[length(curves)+1]]<-data.frame(cohort=co,outcome='All-cause',time=sm$time,group=as.character(sm$strata),risk=1-sm$surv,lower=1-sm$upper,upper=1-sm$lower)
 if(co%in%c('HRS','SHARE')){
  x$ms<-factor(x$status,levels=0:3,labels=c('Censor','Cardiovascular','Non-cardiovascular','Unknown cause'))
  sf<-survfit(Surv(duration,ms)~baseline_heart,data=x,id=person,conf.type='log-log');sm<-summary(sf,censored=TRUE)
  for(j in 2:4)curves[[length(curves)+1]]<-data.frame(cohort=co,outcome=colnames(sm$pstate)[j],time=sm$time,group=as.character(sm$strata),risk=sm$pstate[,j],lower=sm$lower[,j],upper=sm$upper[,j])
 }
}
fwrite(rbindlist(curves),file.path(root,'results/absolute_risk_exact.csv'))
