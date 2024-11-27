#example script to plot violin plots for imaging data from well plates quantified with TissueMaps or Fractal

import glob
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#provide path to metadata file for wells with two columns: well and condition (no column names)
key = pd.read_csv('well_key.txt',sep=' ',header=None,names=['well','condition'])
key = key.set_index('well').to_dict()['condition']

#provide path to csv files quantifying data in each well
infis = glob.glob('csv_dir/*feature-values.csv')
col = '#4af38e'

#load data
for fi in infis:
    df = pd.read_csv(fi,sep=',')
    well = [w for w in key.keys() if w in fi][0]
    df['well'] = well
    df['condition'] = key[well]
    md = pd.read_csv(fi.split('feature-values.csv')[0]+'metadata.csv',sep=',')
    df = df.merge(md,on='mapobject_id')
    if fi == infis[0]:
        all_results = df
    else:
        all_results = pd.concat([all_results,df])

#clean up data
#exclude border nuclei
all_results = all_results[all_results['is_border'] == 0]
#exclude any fields of view that appear out of focus (provide well name, e.g. B02, and x,y coordinates for FOV)
oof_sites = {'wellA':[(0,0)],'wellB':[(0,1)]}
for well in oof_sites:
   for site in oof_sites[well]:
      all_results = all_results[(all_results['well'] != well) | (all_results['well_pos_x'] != oof_sites[well][site[0]]) | (all_results['well_pos_y'] != oof_sites[well][site[1]])]
#all_results = all_results[all_results['nuclei_Count'] == 1] #remove cells with >1 nucleus for analysis if quantifying cellular information (likely mis-segmented)

def plot_violin_intensity(df,cl,seg,stat,stain,coi,control='',norm=False):
    plt.figure(figsize=(3,3))
    if seg != 'nuclei':
        yvar = '{}_{}_Intensity_{}_{}'.format(seg,stat.capitalize(),stat,stain)
    else:
        yvar = 'Intensity_{}_{}'.format(stat,stain)
    if norm:
        df.loc[:,yvar] = (df.loc[:,yvar]-100)/(df.loc[df['condition'] == control,yvar].median()-100) #substract background intensities before calculating L2FC 
    g = sns.violinplot(x='condition',y=yvar,data=df,color=col,order=coi,cut=0)
    xmin, xmax = g.get_xlim()
    if norm:
        plt.hlines(y=1,xmin=xmin, xmax=xmax,colors='lightgrey')
        plt.ylabel('Log2FC {} {}'.format(seg,stain))
        plt.yscale('log',base=2)
    plt.xlabel('')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('figs/{}_violin_{}_{}_{}.pdf'.format(cl,seg,stat,stain),bbox_inches=0,transparent=True)

#filter results & plot
conditions = ['Control','Condition1','Condition2','Condition3'] #conditions to compare (also the order in which to plot conditions)
control = conditions[0]
results = all_results[all_results['condition'].isin(conditions)]
label = 'fig5'
seg = 'nuclei' #name of segmented region
stat = 'mean' #quantified statistic to consider for that segmented region
stain = 'pSR' #name of stain 
plot_violin_intensity(results,label,seg,stat,stain,conditions,control,True)
