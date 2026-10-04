from pathlib import Path
import sys,ast,re,json,hashlib,os
ROOT=Path(os.environ.get('SOURCE_PROJECT_ROOT', '.')).expanduser().resolve();OUT=Path(os.environ.get('ANALYSIS_OUTPUT_ROOT', Path(__file__).resolve().parents[1])).expanduser().resolve()
runtime_root = Path(os.environ.get('RUNTIME_ROOT', ROOT/'output/manuscript_revision_20261004/runtime')).expanduser()
sys.path.insert(0,str(runtime_root))
import numpy as np,pandas as pd,pyreadstat
sys.stdout.reconfigure(encoding='utf-8');(OUT/'private').mkdir(exist_ok=True);(OUT/'results').mkdir(exist_ok=True)
def lit(file,name):
 t=ast.parse(file.read_text(encoding='utf-8-sig'));return next(ast.literal_eval(n.value) for n in t.body if isinstance(n,ast.Assign) and any(isinstance(z,ast.Name) and z.id==name for z in n.targets))
F=lit(ROOT/'合并四库.py','FILES');E=lit(ROOT/'心血管死亡分析.py','EOL');IDS=lit(ROOT/'合并四库.py','ID')
long=pd.read_parquet(ROOT/'分析数据/四库合并长表.parquet');summ=[];endpoints=[];interviews=[];inputs={}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def idstr(s):return s.map(lambda v:str(int(v)) if isinstance(v,(int,float,np.number)) and pd.notna(v) and float(v).is_integer() else str(v).removesuffix('.0'))
def read(path,want):
 _,m=pyreadstat.read_dta(path,metadataonly=True,output_format='dict'); cols=[v for v in m.column_names if want(v.lower())]
 d=pyreadstat.read_dta(path,usecols=cols,output_format='pandas')[0];d.columns=d.columns.str.lower();inputs[path]=sha(path);return d,m
for coh in ['CHARLS','ELSA','HRS','SHARE']:
 g=long[long.cohort.eq(coh)].copy();g['pid']=idstr(g.pid);key=IDS[coh].lower();path=F[coh]
 raw,meta=read(path,lambda v:v in [key,'radyear','radmonth','radsrc','racountry'] or bool(re.match(r'^r\d+(?:iwendm|iwendy|iwm|iwy)$',v)))
 raw['pid']=idstr(raw[key]);assert not raw.pid.duplicated().any()
 dates=[]
 for w in sorted(g.wave.unique()):
  y=f'r{w}iwendy' if coh=='HRS' else f'r{w}iwy';m=f'r{w}iwendm' if coh=='HRS' else f'r{w}iwm'
  d=raw[['pid']].copy();d['wave']=w
  d['year']=raw[y] if y in raw else np.nan;d['month']=raw[m] if m in raw else np.nan;dates.append(d)
 g=g.merge(pd.concat(dates),on=['pid','wave'],validate='one_to_one')
 g['year']=g.year.where(g.year.between(1990,2025));g['month']=g.month.where(g.month.between(1,12))
 g['it']=g.year+(g.month.fillna(6.5)-.5)/12
 # Only true conducted interviews provide covariates and known-alive time.
 g['interview']=g.iwstat.isin([1,2] if coh=='HRS' else [1]) & g.it.notna()
 a=g[g.interview].sort_values(['pid','it','wave']).copy()
 for v in ['hearte','smokev','hibp','diab','rxheart']:a[v]=a[v].where(a[v].isin([0,1]))
 a['age']=a.agey.where(a.agey.between(18,115));a['sex']=a.gender.where(a.gender.isin([1,2]))
 a['edu']=a.educl.where(a.educl.isin([0,1,2,3] if coh=='CHARLS' else [1,2,3]))
 if coh=='CHARLS':a.loc[a.edu.eq(0),'edu']=1
 a['adl']=a.adlfive.where(a.adlfive.between(0,5));a['bmi']=a.bmi.where(a.bmi.between(10,80))
 country=raw.set_index('pid')['racountry'] if 'racountry' in raw else pd.Series(1,index=raw.pid)
 a['country']=a.pid.map(country).fillna(0).astype(int)
 wave_end=a.groupby('wave').it.max().to_dict()
 death_status=g.iwstat.isin([2,3,5,6] if coh=='HRS' else [5,6])
 reports=g[death_status].sort_values(['pid','wave']).drop_duplicates('pid').set_index('pid')
 e,_=read(E[coh][0],lambda v:v in [key,'raxyear','raxmonth','raxseason','ragcod']);e['pid']=idstr(e[key]);assert not e.pid.duplicated().any();e=e.set_index('pid')
 r=raw.set_index('pid');ss={'cohort':coh,'source_n':int(g.pid.nunique()),'actual_interview_n':int(a.pid.nunique())}
 # Retain all cleaned interview rows; eligibility is locked below at age >=50.
 interviews.append(a[['cohort','pid','wave','it','year','month','age','sex','edu','smokev','hibp','adl','bmi','hearte','rxheart','diab','dep_bin','country']])
 person=[]
 for pid,h in a.groupby('pid',sort=False):
  first=h.iloc[0];last=float(h.it.max());year=r.at[pid,'radyear'] if 'radyear' in r else np.nan
  month=r.at[pid,'radmonth'] if 'radmonth' in r else np.nan;source='core_death_date'
  if not pd.notna(year) or not 1990<=year<=2025:
   year=e.at[pid,'raxyear'] if pid in e.index and 'raxyear' in e else np.nan
   month=e.at[pid,'raxmonth'] if pid in e.index and 'raxmonth' in e else np.nan;source='eol_death_date'
  known_year=pd.notna(year) and 1990<=year<=2025;known_month=pd.notna(month) and 1<=month<=12
  died=known_year or pid in reports.index or pid in e.index
  cause=e.at[pid,'ragcod'] if pid in e.index else np.nan
  status=1 if cause==2 else 2 if cause in [1,3] else 3
  lo=hi=np.nan;conflict=False
  if died:
   if known_year:
    lo=float(year)+(float(month)-1)/12 if known_month else float(year)
    hi=lo+1/12 if known_month else float(year)+1
   else:
    source='death_notification_interval';lo=last
    if pid in reports.index:
     w=reports.at[pid,'wave'];hi=wave_end.get(w,np.nan)
     if pd.notna(hi):hi=float(hi)
   if pd.notna(lo) and pd.notna(hi):
    conflict=last>=hi
    lo=max(lo,last) # Same-month deaths occur after the completed interview.
   if pd.isna(hi) or hi<=lo:conflict=True
  person.append({'cohort':coh,'pid':pid,'died':int(died),'cause':status if died else 0,'death_lo':lo,'death_hi':hi,'last_alive':last,'death_source':source if died else 'last_interview','date_conflict':int(conflict),'first_interview':float(first.it)})
 p=pd.DataFrame(person);p.to_parquet(OUT/'private'/f'endpoints_{coh}.parquet',index=False)
 ss.update(death_notified=int(p.died.sum()),invalid_death_interval=int(p.date_conflict.sum()),dated_deaths=int((p.died.eq(1)&p.death_source.ne('death_notification_interval')).sum()))
 endpoints.append(p);summ.append(ss);print(ss,flush=True)
iv=pd.concat(interviews,ignore_index=True);ep=pd.concat(endpoints,ignore_index=True)
iv.to_parquet(OUT/'private/interviews.parquet',index=False);ep.to_parquet(OUT/'private/endpoints.parquet',index=False)
inputs[str(ROOT/'分析数据/四库合并长表.parquet')]=sha(ROOT/'分析数据/四库合并长表.parquet')
(OUT/'input_hashes.json').write_text(json.dumps(inputs,ensure_ascii=False,indent=2),encoding='utf-8')
pd.DataFrame(summ).to_csv(OUT/'results/source_endpoints_summary.csv',index=False)
# Baseline: first actual age-eligible interview, then exclusions; do not delay baseline until disease is observed.
base=iv[iv.age.ge(50)].sort_values(['cohort','pid','it','wave']).drop_duplicates(['cohort','pid']).copy()
base=base.merge(ep,on=['cohort','pid'],validate='one_to_one');base['complete']=base[['sex','edu','smokev','hibp','adl']].notna().all(axis=1)
base.to_parquet(OUT/'private/baseline_all.parquet',index=False)
ss=[]
for c,g in base.groupby('cohort'):
 ss.append(dict(cohort=c,age50_n=len(g),exposure_missing=int(g.hearte.isna().sum()),invalid_endpoint=int((g.hearte.notna()&g.date_conflict.eq(1)).sum()),complete_covariate=int(g.complete.sum()),bmi_missing=int(g.bmi.isna().sum()),baseline_min=float(g.it.min()),baseline_max=float(g.it.max()),age_min=float(g.age.min()),age_max=float(g.age.max())))
pd.DataFrame(ss).to_csv(OUT/'results/baseline_initial_audit.csv',index=False);print(pd.DataFrame(ss).to_string(index=False))
