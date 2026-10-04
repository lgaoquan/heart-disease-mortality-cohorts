from pathlib import Path
import sys,json,hashlib,os
OUT=Path(os.environ.get('ANALYSIS_OUTPUT_ROOT', Path(__file__).resolve().parents[1])).expanduser().resolve();ROOT=Path(os.environ.get('SOURCE_PROJECT_ROOT', OUT.parent.parent)).expanduser().resolve()
runtime_root = Path(os.environ.get('RUNTIME_ROOT', ROOT/'output/manuscript_revision_20261004/runtime')).expanduser()
sys.path.insert(0,str(runtime_root))
import pandas as pd,numpy as np
sys.stdout.reconfigure(encoding='utf-8')
audit_root = Path(os.environ.get('AUDIT_ROOT', ROOT/'output/manuscript_revision_20261004/code_audit')).expanduser()
sys.path.insert(0,str(audit_root))
from validated_primitives import self_test
self_test()
b=pd.read_csv(OUT/'private/baseline.csv');mid=pd.read_csv(OUT/'private/intervals_mid_5y.csv')
b['duration']=b.duration.round(8)
checks={};checks['event_person_counts_match']=bool(mid.groupby('cohort').event.sum().equals(b.groupby('cohort').event.sum()))
checks['no_repeated_deaths']=bool(mid.groupby('person').event.sum().le(1).all())
checks['starts_at_zero']=bool(mid.groupby('person').start.min().eq(0).all())
checks['no_nonpositive_intervals']=bool(mid.stop.gt(mid.start).all())
checks['no_interval_beyond_5y']=bool(mid.stop.le(5+1e-8).all())
last=mid.sort_values(['person','start']).drop_duplicates('person',keep='last')
checks['all_deaths_terminal']=int(last.event.sum())==int(mid.event.sum())
r=pd.read_csv(OUT/'results/model_results.csv');mc=pd.read_csv(OUT/'results/primary_counts.csv')
checks['primary_model_counts_agree']=all(int(r[(r.cohort==z.cohort)&(r.analysis=='primary')&(r.outcome=='all')].iloc[0]['n'])==z.n and int(r[(r.cohort==z.cohort)&(r.analysis=='primary')&(r.outcome=='all')].iloc[0]['events'])==z.deaths for z in mc.itertuples())
checks['all_primary_estimates_finite']=bool(np.isfinite(r[r.analysis.eq('primary')][['hr','lower','upper','beta','se']]).all().all())
checks['all_lower_hr_upper']=bool(((r.lower<r.hr)&(r.hr<r.upper)).all())
# Independently calculate unique-person all-cause and disjoint-cause risk at 5 y.
rr=pd.read_csv(OUT/'results/absolute_risk_5y.csv');comparisons=[]
for (co,ex),g in b[b.complete].groupby(['cohort','baseline_heart']):
 s=1.;cif=np.zeros(3)
 counts=g.groupby('duration').status.agg(list)
 y=len(g)
 for t,st in counts.items():
  ev=np.array([st.count(j) for j in [1,2,3]])
  cif+=s*ev/y;s*=1-sum(ev)/y;y-=len(st)
  if abs(s+cif.sum()-1)>1e-9:raise AssertionError('AJ probability conservation failed')
 q=rr[(rr.cohort==co)&rr['group'].eq('baseline_heart='+str(int(ex)))]
 diff=abs(q[q.outcome.eq('All-cause')].iloc[0].risk-(1-s))
 comparisons.append({'cohort':co,'exposure':int(ex),'outcome':'All-cause','absolute_difference':diff})
 if co in ['HRS','SHARE']:
  for j,name in enumerate(['Cardiovascular','Non-cardiovascular','Unknown cause']):
   diff=abs(q[q.outcome.eq(name)].iloc[0].risk-cif[j]);comparisons.append({'cohort':co,'exposure':int(ex),'outcome':name,'absolute_difference':diff})
checks['independent_AJ_and_KM_agree']=bool(max(x['absolute_difference'] for x in comparisons)<1e-8)
checks['original_files_unchanged']=True
old_path = Path(os.environ.get('INPUT_HASHES', ROOT/'output/manuscript_revision_20261004/evidence/input_hashes.json')).expanduser()
old=json.loads(old_path.read_text(encoding='utf-8'))
for p,h in old.items():
 if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h:checks['original_files_unchanged']=False
record={'checks':checks,'all_pass':all(checks.values()),'independent_risk_comparison':comparisons,'limitations':['Measurement, sampling, timing assumptions and residual confounding remain substantive limitations.','Synthetic checks and numerical agreement do not establish ethics or author approval.']}
(OUT/'results/validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record,indent=2));assert all(checks.values())
