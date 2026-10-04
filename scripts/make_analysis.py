from pathlib import Path
import os
import sys,json
OUT=Path(os.environ.get('ANALYSIS_OUTPUT_ROOT', Path(__file__).resolve().parents[1])).expanduser().resolve();ROOT=Path(os.environ.get('SOURCE_PROJECT_ROOT', OUT.parent.parent)).expanduser().resolve()
sys.path.insert(0,str(ROOT/'output/manuscript_revision_20261004/runtime'))
import pandas as pd,numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base=pd.read_parquet(OUT/'private/baseline_all.parquet');iv=pd.read_parquet(OUT/'private/interviews.parquet')
base['person']=np.arange(len(base))+1
flows=[];bi=[];missing=[];diagnostics=[]
for scenario in ['mid','early','late']:
 b=base.copy();b['terminal']=b.last_alive
 d=b.died.eq(1)&b.date_conflict.eq(0)
 when=(b.death_lo+b.death_hi)/2 if scenario=='mid' else b.death_lo+np.minimum(1/365.25,(b.death_hi-b.death_lo)/2) if scenario=='early' else b.death_hi-1e-8
 b.loc[d,'terminal']=when[d]
 for c,g in b.groupby('cohort'):
  if scenario!='mid':continue
  masks=[g.hearte.notna(),g.date_conflict.eq(0),g.terminal.gt(g.it)]
  sel=np.ones(len(g),dtype=bool);counts={'cohort':c,'age_eligible':len(g)}
  for name,mask in zip(['missing_exposure','invalid_death_interval','no_positive_followup'],masks):counts[name]=int((sel&~mask).sum());sel&=mask
  counts['eligible']=int(sel.sum());counts['missing_covariates']=int((sel&~g.complete).sum());counts['primary_n']=int((sel&g.complete).sum());flows.append(counts)
 b=b[b.hearte.notna()&b.date_conflict.eq(0)&b.terminal.gt(b.it)].copy()
 b['base_time']=b.it;b['baseline_heart']=b.hearte;b['baseline_med']=np.where(b.hearte.eq(0),0,np.where(b.rxheart.eq(0),1,np.where(b.rxheart.eq(1),2,np.nan)))
 for horizon in ([5,3] if scenario=='mid' else [5]):
  h=b.copy();h['end']=np.minimum(h.terminal,h.it+horizon);h['duration']=h.end-h.it
  h['status']=np.where(h.died.eq(1)&h.terminal.le(h.it+horizon),h.cause,0).astype(int)
  h['event']=(h.status>0).astype(int)
  h['date_known']=h.death_source.ne('death_notification_interval').astype(int)
  # All-cause exactly one event per person; future deaths become administrative censoring.
  kept=['cohort','pid','person','base_time','end','duration','status','event','age','sex','edu','smokev','hibp','adl','bmi','baseline_heart','baseline_med','country','complete','date_known','death_source','diab','dep_bin']
  h=h[kept].copy();h['base_year']=np.floor(h.base_time)
  if scenario=='mid' and horizon==5:
   h.drop(columns=['pid']).to_csv(OUT/'private/baseline.csv',index=False)
   h.to_parquet(OUT/'private/final_baseline.parquet',index=False)
   for c,g in h.groupby('cohort'):
    for expo,x in g.groupby('baseline_heart'):
     for v in ['age','sex','edu','smokev','hibp','adl','bmi','diab','baseline_med']:
      missing.append({'cohort':c,'baseline_heart':int(expo),'variable':v,'denominator':len(x),'missing_n':int(x[v].isna().sum()),'missing_percent':100*x[v].isna().mean()})
   for c,g in h[h.complete].groupby('cohort'):
    diagnostics.append({'cohort':c,'n':len(g),'deaths':int(g.event.sum()),'cv':int(g.status.eq(1).sum()),'noncv':int(g.status.eq(2).sum()),'unknown':int(g.status.eq(3).sum()),'person_years':float(g.duration.sum()),'median_followup':float(g.duration.median()),'approx_interval_deaths':int((g.event.eq(1)&g.death_source.eq('death_notification_interval')).sum()),'baseline_min':float(g.base_time.min()),'baseline_max':float(g.base_time.max()),'exposed':int(g.baseline_heart.sum())})
  # Prepare time-varying risk intervals; baseline covariates remain fixed.
  obs=iv[['cohort','pid','wave','it','hearte','rxheart']].merge(h,on=['cohort','pid'],validate='many_to_one')
  obs=obs[obs.it.ge(obs.base_time)&obs.it.lt(obs.end)].sort_values(['person','it','wave']).drop_duplicates(['person','it'],keep='last')
  obs['heart']=obs.groupby('person',sort=False).hearte.ffill()
  obs['rx']=obs.groupby('person',sort=False).rxheart.ffill()
  obs['start']=obs.it-obs.base_time;obs['stop']=obs.groupby('person',sort=False).it.shift(-1).fillna(obs.end)-obs.base_time
  obs['stop']=np.minimum(obs.stop,obs.duration)
  obs['event']=obs.event*obs.stop.ge(obs.duration-1e-8);obs['status']=obs.status*obs.stop.ge(obs.duration-1e-8)
  obs['late_period']=(obs.start>=2.5).astype(int)
  assert obs.heart.notna().all() and obs.stop.gt(obs.start).all()
  assert obs.groupby('person').event.sum().le(1).all()
  assert obs.groupby('person').start.min().eq(0).all()
  assert obs.groupby('person').stop.max().sort_index().equals(h.set_index('person').duration.sort_index())
  assert int(obs.event.sum())==int(h.event.sum()) and obs.person.nunique()==len(h)
  outcols=['cohort','person','start','stop','status','event','heart','baseline_heart','baseline_med','age','sex','edu','smokev','hibp','adl','bmi','country','base_year','complete','date_known']
  name=f'{scenario}_{horizon}y';obs[outcols].to_csv(OUT/'private'/f'intervals_{name}.csv',index=False)
  print(name,len(h),'persons',len(obs),'intervals',int(obs.event.sum()),'unique deaths',flush=True)
pd.DataFrame(flows).to_csv(OUT/'results/flow.csv',index=False)
pd.DataFrame(missing).to_csv(OUT/'results/missingness.csv',index=False)
pd.DataFrame(diagnostics).to_csv(OUT/'results/primary_counts.csv',index=False)
(OUT/'results/interval_validation.json').write_text(json.dumps({'unique_death_per_person':True,'first_interval_retained':True,'no_risk_after_endpoint':True,'administrative_censoring_resets_events':True,'identical_participant_and_interval_events':True,'construction':'independent revision; all original files retained'},indent=2),encoding='utf-8')
