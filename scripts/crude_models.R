suppressPackageStartupMessages({library(survival);library(data.table)})
root<-Sys.getenv('ANALYSIS_OUTPUT_ROOT', unset='analysis_output');d<-fread(file.path(root,'private/intervals_mid_5y.csv'));r<-fread(file.path(root,'results/model_results.csv'));r<-r[analysis!='unadjusted'];new<-list()
for(co in c('CHARLS','ELSA','HRS','SHARE')){
 x<-d[cohort==co & complete==TRUE];x[,c('start','stop'):=.(round(start,8),round(stop,8))]
 f<-coxph(Surv(start,stop,event)~heart,data=x,cluster=person,robust=TRUE,ties='efron');bt<-unname(coef(f)['heart']);se<-sqrt(vcov(f)['heart','heart'])
 new[[length(new)+1]]<-data.frame(cohort=co,analysis='unadjusted',outcome='all',term='heart',n=uniqueN(x$person),intervals=nrow(x),events=sum(x$event),beta=bt,se=se,hr=exp(bt),lower=exp(bt-1.96*se),upper=exp(bt+1.96*se),p=2*pnorm(-abs(bt/se)),warnings='')
}
fwrite(rbindlist(c(list(r),new),fill=TRUE),file.path(root,'results/model_results.csv'))
