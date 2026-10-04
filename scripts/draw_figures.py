from pathlib import Path
import os
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Polygon
from matplotlib.ticker import PercentFormatter,ScalarFormatter
from PIL import Image
O=Path(os.environ.get('ANALYSIS_OUTPUT_ROOT', Path(__file__).resolve().parents[1])).expanduser().resolve();R=O/'results';F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
BLUE='#325A70';GREY='#7C8991';DARK='#222D35';CO=['CHARLS','ELSA','HRS','SHARE']
flow=pd.read_csv(R/'flow.csv').set_index('cohort');counts=pd.read_csv(R/'primary_counts.csv').set_index('cohort');src=pd.read_csv(R/'source_endpoints_summary.csv').set_index('cohort')
models=pd.read_csv(R/'model_results.csv');meta=pd.read_csv(R/'meta_results.csv');curves=pd.read_csv(R/'absolute_risk_curves.csv');exact=pd.read_csv(R/'absolute_risk_exact.csv')
def save(fig,name):
 fig.savefig(F/(name+'.pdf'))
 fig.savefig(F/(name+'.png'),dpi=600)
 im=Image.open(F/(name+'.png')).convert('RGB');im.save(F/(name+'.png'),dpi=(600,600));im.save(F/(name+'.tiff'),dpi=(600,600),compression='tiff_lzw')
 im.thumbnail((1400,1800));im.save(F/(name+'_review.png'));plt.close(fig)
# Flow diagram with an aligned numerical structure and explicitly sequential exclusions.
fig=plt.figure(figsize=(170/25.4,155/25.4));ax=fig.add_axes([.02,.02,.96,.95]);ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
xx=np.linspace(.395,.915,4)
for x,c in zip(xx,CO):ax.text(x,.96,c,ha='center',weight='bold',fontsize=9,color=BLUE)
rows=[('Source participants',lambda c:int(src.loc[c,'source_n']),False),('First interview at age ≥50',lambda c:int(flow.loc[c,'age_eligible']),False),('Excluded: exposure missing',lambda c:int(flow.loc[c,'missing_exposure']),True),('Excluded: inconsistent death dates',lambda c:int(flow.loc[c,'invalid_death_interval']),True),('Excluded: no positive follow-up',lambda c:int(flow.loc[c,'no_positive_followup']),True),('Eligible with follow-up',lambda c:int(flow.loc[c,'eligible']),False),('Excluded: adjustment data missing',lambda c:int(flow.loc[c,'missing_covariates']),True),('Primary analysis',lambda c:int(flow.loc[c,'primary_n']),False),('Unique deaths within 5 years',lambda c:int(counts.loc[c,'deaths']),False)]
for i,(label,fun,excluded) in enumerate(rows):
 y=.88-i*.087
 ax.text(.015,y,label,ha='left',va='center',fontsize=7.1 if excluded else 8,weight='normal' if excluded else 'bold',color=GREY if excluded else DARK)
 for x,c in zip(xx,CO):
  if not excluded:ax.add_patch(Rectangle((x-.073,y-.027),.146,.054,facecolor='#F2F5F6',edgecolor='#C3CCD1',lw=.55))
  ax.text(x,y,str(fun(c)),ha='center',va='center',fontsize=8.4,color=GREY if excluded else DARK,weight='normal' if excluded else 'bold')
 ax.plot([.01,.995],[y-.044,y-.044],color='#E0E5E8',lw=.45)
ax.text(.015,.04,'Exclusions are sequential. Adjustment data: sex, education, ever smoking, hypertension and ADL.\nAge is required for eligibility; follow-up starts at the first age-eligible interview.',fontsize=7,color=DARK,va='bottom')
save(fig,'Figure1_flow')
# Forest plot: one coherent numerical column, independent panel axes documented.
fig=plt.figure(figsize=(170/25.4,180/25.4));gs=fig.add_gridspec(3,3,width_ratios=[1.25,2.3,1.25],height_ratios=[1.45,1,1],left=.035,right=.98,top=.95,bottom=.11,hspace=.52,wspace=.04)
for panel,(outcome,title) in enumerate([('all','A  All-cause mortality'),('cv','B  Cardiovascular mortality'),('noncv','C  Non-cardiovascular mortality')]):
 q=models[models.analysis.eq('primary')&models.outcome.eq(outcome)&models.cohort.isin(CO if outcome=='all' else ['HRS','SHARE'])].set_index('cohort')
 cohort=[c for c in CO if c in q.index];ys=np.arange(len(cohort),0,-1)
 la=fig.add_subplot(gs[panel,0]);a=fig.add_subplot(gs[panel,1]);ra=fig.add_subplot(gs[panel,2]);la.axis('off');ra.axis('off')
 ylow=-1 if outcome=='all' else .3;yhigh=len(cohort)+.85
 for z in [la,a,ra]:z.set_ylim(ylow,yhigh)
 la.text(0,yhigh,title,fontsize=9,weight='bold',va='bottom',transform=la.get_yaxis_transform(),clip_on=False)
 la.text(0,len(cohort)+.42,'Cohort (deaths)',weight='bold',fontsize=7.7);ra.text(.02,len(cohort)+.42,'HR (95% CI)',weight='bold',fontsize=7.7)
 for y,c in zip(ys,cohort):
  row=q.loc[c];la.text(0,y,f'{c} ({int(row.events)})',va='center',fontsize=8)
  a.errorbar(row.hr,y,xerr=[[row.hr-row.lower],[row.upper-row.hr]],fmt='s',ms=4,color=BLUE,lw=.9,capsize=2)
  ra.text(.02,y,f'{row.hr:.2f} ({row.lower:.2f}–{row.upper:.2f})',va='center',fontsize=8)
 if outcome=='all':
  m=meta[meta.analysis.eq('primary')&meta.outcome.eq('all')].iloc[0];y=0
  la.text(0,y,'Random effects',va='center',fontsize=8,weight='bold');a.add_patch(Polygon([[m.lower,y],[m.hr,y+.14],[m.upper,y],[m.hr,y-.14]],color=BLUE));ra.text(.02,y,f'{m.hr:.2f} ({m.lower:.2f}–{m.upper:.2f})',va='center',fontsize=8,weight='bold')
  a.plot([m.pi_lower,m.pi_upper],[-.55,-.55],color=GREY,lw=1,ls='--');la.text(0,-.55,'Prediction interval',va='center',fontsize=7,color=GREY);ra.text(.02,-.55,f'{m.pi_lower:.2f}–{m.pi_upper:.2f}',va='center',fontsize=7,color=GREY)
 a.set_xscale('log');a.set_xlim((.5,4.2) if outcome=='all' else (.6,4.2));a.set_xticks([.5,1,2,4] if outcome=='all' else [1,2,4]);a.xaxis.set_major_formatter(ScalarFormatter());a.minorticks_off();a.set_yticks([]);a.spines['left'].set_visible(False);a.axvline(1,color=GREY,lw=.7,ls=':');a.set_xlabel('Hazard ratio (log scale)',fontsize=7.5)
fig.text(.035,.012,'Panel-specific axes. Cause-specific analyses are restricted to HRS and SHARE; no pooled cause-specific diamond is shown.',fontsize=6.5,color=DARK)
save(fig,'Figure2_forest')
def risk_panels(pairs,name,height):
 fig,axs=plt.subplots(2,2,figsize=(170/25.4,height/25.4));fig.subplots_adjust(left=.11,right=.97,top=.91,bottom=.20,wspace=.26,hspace=.92)
 for j,(ax,(co,outcome)) in enumerate(zip(axs.flat,pairs)):
  g=curves[curves.cohort.eq(co)&curves.outcome.eq(outcome)]
  for expo,color,lab in [(0,GREY,'No heart disease'),(1,BLUE,'Heart disease')]:
   q=g[g['group'].eq(f'baseline_heart={expo}')].sort_values('time');z=exact[exact.cohort.eq(co)&exact.outcome.eq(outcome)&exact['group'].eq(f'baseline_heart={expo}')].sort_values('time')
   ax.step(np.r_[0,z.time],np.r_[0,z.risk],where='post',color=color,lw=1.1,label=lab);ax.fill_between(np.r_[0,z.time],np.r_[0,z.lower],np.r_[0,z.upper],step='post',alpha=.14,color=color,lw=0)
   for ti,x in enumerate([0,1,3,5]):
    z=q[np.isclose(q.time,x)].iloc[0];ax.text(x/5,-.40-expo*.14,str(int(z.at_risk)),transform=ax.transAxes,ha='center',fontsize=6.8,color=color)
  ax.text(-.03,-.27,'Number at risk',transform=ax.transAxes,fontsize=6.7,ha='left')
  ax.set_title(chr(65+j)+'  '+co+(' | '+outcome if outcome!='All-cause' else ''),loc='left',fontsize=8.5,weight='bold',pad=7)
  ax.set_ylim(0,.35 if outcome=='All-cause' else .16);ax.set_yticks(np.arange(0,.351,.05) if outcome=='All-cause' else [0,.05,.10,.15]);ax.set_xlim(0,5);ax.set_xticks([0,1,3,5]);ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.set_xlabel('Years since baseline',fontsize=7.5);ax.set_ylabel('Cumulative mortality' if outcome=='All-cause' else 'Cumulative incidence',fontsize=7.5);ax.grid(axis='y',color='#E0E5E8',lw=.4)
 handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=2,frameon=False,fontsize=8,bbox_to_anchor=(.55,.988))
 save(fig,name)
risk_panels([(c,'All-cause') for c in CO],'Figure3_cumulative_mortality',175)
risk_panels([(c,o) for c in ['HRS','SHARE'] for o in ['Cardiovascular','Non-cardiovascular']],'FigureS1_cause_specific_risks',175)
print('Exported four checked-layout figures in PDF, PNG and RGB LZW TIFF.')
