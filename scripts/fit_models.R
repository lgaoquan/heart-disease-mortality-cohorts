suppressPackageStartupMessages({library(survival);library(data.table);library(splines);library(metafor);library(jsonlite)})
root <- Sys.getenv('ANALYSIS_OUTPUT_ROOT', unset='analysis_output')
out <- file.path(root,'results');dir.create(out,showWarnings=FALSE)
d <- fread(file.path(root,'private/intervals_mid_5y.csv'))
b <- fread(file.path(root,'private/baseline.csv'))
d[,c('start','stop'):=.(round(start,8),round(stop,8))];b[,duration:=round(duration,8)]
results<-list(); failures<-list();phs<-list(); interaction<-list();meds<-list(); idx<-0L
capture_fit <- function(g,analysis,outcome='all',exposure='heart',adjust='full',diagnostic=FALSE) {
 g<-copy(g);g[,c('start','stop'):=.(round(start,8),round(stop,8))];g<-g[stop-start>1e-7];co<-as.character(g$cohort[1]);g[,ev:=if(outcome=='all') event else if(outcome=='cv') as.integer(status==1) else if(outcome=='noncv') as.integer(status==2) else if(outcome=='cv_unknown') as.integer(status%in%c(1,3)) else as.integer(status%in%c(2,3))]
 if(length(unique(g$person))<100 || sum(g$ev)<20) {failures[[length(failures)+1]]<<-data.frame(cohort=co,analysis=analysis,outcome=outcome,reason='Fewer than 100 people or 20 events');return(NULL)}
 g[,sex:=factor(sex)];g[,edu:=factor(edu)];g[,smokev:=factor(smokev)];g[,hibp:=factor(hibp)]
 rhs<-paste(exposure,'+ ns(age,df=3) + sex + base_year')
 if(adjust%in%c('full','bmi','missing','timeflex'))rhs<-paste(rhs,'+ edu + smokev + hibp + adl')
 if(adjust=='bmi')rhs<-paste(rhs,'+ ns(bmi,df=3)')
 if(adjust=='missing')rhs<-paste(rhs,'+ adl_missing')
 if(adjust=='timeflex')rhs<-paste(rhs,'+ period_late:(ns(age,df=3) + sex + edu + smokev + hibp + adl + base_year)')
 if(co=='SHARE')rhs<-paste(rhs,'+ strata(country)')
 fm<-as.formula(paste('Surv(start,stop,ev) ~',rhs))
 warnings<-character()
 fit<-tryCatch(withCallingHandlers(coxph(fm,data=g,cluster=person,robust=TRUE,ties='efron',x=diagnostic,model=FALSE,control=coxph.control(iter.max=50)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error')) {failures[[length(failures)+1]]<<-data.frame(cohort=co,analysis=analysis,outcome=outcome,reason=conditionMessage(fit));return(NULL)}
 if(length(warnings)>0)failures[[length(failures)+1]]<<-data.frame(cohort=co,analysis=analysis,outcome=outcome,reason=paste(warnings,collapse='; '))
 names_keep<-if(exposure=='factor(baseline_med)')grep('baseline_med',names(coef(fit)),value=TRUE) else if(exposure=='heart + heart_late')c('heart','heart_late') else exposure
 for(n in names_keep){
  z<-coef(fit)[n];se<-sqrt(diag(vcov(fit)))[n]
  if(is.na(z)||!is.finite(se))next
  idx<<-idx+1L;results[[idx]]<<-data.frame(cohort=co,analysis=analysis,outcome=outcome,term=n,n=uniqueN(g$person),intervals=nrow(g),events=sum(g$ev),beta=unname(z),se=unname(se),hr=exp(z),lower=exp(z-1.96*se),upper=exp(z+1.96*se),p=2*pnorm(-abs(z/se)),warnings=paste(warnings,collapse='; '),row.names=NULL)
 }
 if(diagnostic){
  zz<-tryCatch(cox.zph(fit,transform='km'),error=function(e)e)
  if(!inherits(zz,'error')){pp<-as.data.frame(zz$table);pp$term<-rownames(pp);pp$cohort<-co;pp$analysis<-analysis;pp$outcome<-outcome;phs[[length(phs)+1]]<<-pp}
 }
 cat(co,analysis,outcome,'N',uniqueN(g$person),'E',sum(g$ev),'\n');flush.console()
 return(fit)
}
# Main sample is complete for common baseline adjustment variables, without future BMI measurements.
for(co in c('CHARLS','ELSA','HRS','SHARE')){
 x<-d[cohort==co & complete==TRUE]
 for(oc in c('all','cv','noncv'))capture_fit(x,'primary',oc,diagnostic=TRUE)
 capture_fit(d[cohort==co & !is.na(sex)],'age_sex',adjust='minimal')
 capture_fit(x,'baseline_exposure',exposure='baseline_heart')
 capture_fit(x[age>=60],'age60')
 capture_fit(x[date_known==1],'known_date_only')
 if(co!='ELSA')capture_fit(x[!is.na(bmi)],'bmi_extension',adjust='bmi')
 fs<-capture_fit(x,'sex_interaction',exposure='heart * sex')
 if(!is.null(fs)){
  cv<-vcov(fs);cf<-coef(fs);tn<-'heart:sex2';if(tn%in%names(cf)){
   mu<-c(cf['heart'],cf['heart']+cf[tn]);ss<-c(sqrt(cv['heart','heart']),sqrt(cv['heart','heart']+cv[tn,tn]+2*cv['heart',tn]))
   meds[[length(meds)+1]]<-data.frame(cohort=co,sex=c('Men','Women'),hr=exp(mu),lower=exp(mu-1.96*ss),upper=exp(mu+1.96*ss),interaction_p=2*pnorm(-abs(cf[tn]/sqrt(cv[tn,tn]))))
  }
 }
 # Test and estimate an early/late heart-disease association under non-proportional hazards.
 gx<-as.data.frame(x);gx$ev<-gx$event
 sp<-survSplit(Surv(start,stop,ev)~.,data=gx,cut=2.5,episode='period')
 sp$event<-sp$ev;sp$status<-sp$status*sp$ev
 stopifnot(sum(sp$event)==sum(x$event),all(tapply(sp$event,sp$person,sum)<=1))
 sp$heart_late<-sp$heart*(sp$start>=2.5)
 sp$period_late<-as.integer(sp$start>=2.5)
 capture_fit(as.data.table(sp),'nuisance_timeflex',adjust='timeflex')
 f<-capture_fit(as.data.table(sp),'time_interaction',exposure='heart + heart_late')
 if(!is.null(f)){
  vv<-vcov(f);cc<-coef(f);v<-c(early=cc['heart'],late=cc['heart']+cc['heart_late']);se<-c(sqrt(vv['heart','heart']),sqrt(vv['heart','heart']+vv['heart_late','heart_late']+2*vv['heart','heart_late']))
  interaction[[length(interaction)+1]]<-data.frame(cohort=co,period=c('0-2.5 years','2.5-5 years'),hr=exp(v),lower=exp(v-1.96*se),upper=exp(v+1.96*se),interaction_p=2*pnorm(-abs(cc['heart_late']/sqrt(vv['heart_late','heart_late']))))
 }
 xm<-x[!is.na(baseline_med)];if(length(unique(xm$baseline_med))==3)capture_fit(xm,'baseline_medication',exposure='factor(baseline_med)')
 # Missing-category sensitivity is explicitly separate from complete-case primary inference.
 xm<-copy(d[cohort==co & !is.na(sex)]);for(v in c('edu','smokev','hibp'))set(xm,which(is.na(xm[[v]])),v,9)
 xm[,adl_missing:=as.integer(is.na(adl))];xm[is.na(adl),adl:=0]
 capture_fit(xm,'missing_category',adjust='missing')
 if(co%in%c('HRS','SHARE'))for(oc in c('cv_unknown','noncv_unknown'))capture_fit(x,'unknown_reassigned',oc)
}
for(sc in c('mid_3y','early_5y','late_5y')){
 z<-fread(file.path(root,paste0('private/intervals_',sc,'.csv')))
 for(co in c('CHARLS','ELSA','HRS','SHARE'))capture_fit(z[cohort==co & complete==TRUE],sc)
}
r<-rbindlist(results,fill=TRUE);fwrite(r,file.path(out,'model_results.csv'))
fwrite(rbindlist(failures,fill=TRUE),file.path(out,'model_warnings.csv'))
fwrite(rbindlist(phs,fill=TRUE),file.path(out,'proportional_hazards.csv'))
fwrite(rbindlist(interaction,fill=TRUE),file.path(out,'time_interaction.csv'))
fwrite(rbindlist(meds,fill=TRUE),file.path(out,'sex_interaction.csv'))
# REML with modified Hartung-Knapp SE, t intervals; no four-cohort cause-specific pooling.
metas<-list()
for(an in unique(r$analysis))for(oc in unique(r$outcome)){
 q<-r[analysis==an & outcome==oc & term%in%c('heart','baseline_heart') & warnings=='']
 if(oc!='all')q<-q[cohort%in%c('HRS','SHARE')]
 k<-nrow(q);if(k<2)next
 m<-tryCatch(rma.uni(yi=q$beta,sei=q$se,method='REML',test='knha'),error=function(e)NULL);if(is.null(m))next
 w<-1/(q$se^2+m$tau2);mu<-sum(w*q$beta)/sum(w);scale<-max(1,sum(w*(q$beta-mu)^2)/(k-1));se<-sqrt(scale/sum(w));cr<-qt(.975,k-1)
 pi<-if(k>2)qt(.975,k-2)*sqrt(m$tau2+se^2) else NA_real_
 metas[[length(metas)+1]]<-data.frame(analysis=an,outcome=oc,k=k,n=sum(q$n),events=sum(q$events),beta=mu,se=se,hr=exp(mu),lower=exp(mu-cr*se),upper=exp(mu+cr*se),I2=m$I2,tau2=m$tau2,Q=m$QE,Q_p=m$QEp,pi_lower=exp(mu-pi),pi_upper=exp(mu+pi))
}
fwrite(rbindlist(metas),file.path(out,'meta_results.csv'))
# Unique-person absolute risks, baseline exposure. No time-varying exposure grouping.
curves<-list();risk5<-list()
for(co in c('CHARLS','ELSA','HRS','SHARE')){
 x<-as.data.frame(b[cohort==co & complete==TRUE]);x$baseline_heart<-factor(x$baseline_heart,levels=c(0,1))
 sf<-survfit(Surv(duration,event)~baseline_heart,data=x,conf.type='log-log')
 sm<-summary(sf,times=seq(0,5,.05),extend=FALSE)
 cc<-data.frame(cohort=co,outcome='All-cause',time=sm$time,group=as.character(sm$strata),at_risk=sm$n.risk,risk=1-sm$surv,lower=1-sm$upper,upper=1-sm$lower)
 curves[[length(curves)+1]]<-cc;risk5[[length(risk5)+1]]<-cc[abs(cc$time-5)<1e-8,]
 if(co%in%c('HRS','SHARE')){
  x$ms<-factor(x$status,levels=0:3,labels=c('Censor','Cardiovascular','Non-cardiovascular','Unknown cause'))
  sf<-survfit(Surv(duration,ms)~baseline_heart,data=x,id=person,conf.type='log-log')
  sm<-summary(sf,times=seq(0,5,.05),extend=FALSE)
  for(j in 2:4){cc<-data.frame(cohort=co,outcome=colnames(sm$pstate)[j],time=sm$time,group=as.character(sm$strata),at_risk=sm$n.risk[,1],risk=sm$pstate[,j],lower=sm$lower[,j],upper=sm$upper[,j]);curves[[length(curves)+1]]<-cc;risk5[[length(risk5)+1]]<-cc[abs(cc$time-5)<1e-8,]}
 }
}
fwrite(rbindlist(curves),file.path(out,'absolute_risk_curves.csv'));fwrite(rbindlist(risk5),file.path(out,'absolute_risk_5y.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
cat('DONE: models, diagnostics, REML, unique-person risks exported.\n')
