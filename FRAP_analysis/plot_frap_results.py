import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_curves(frap_results,conditions,colours,maxT,figname,outdir):
   plt.figure(figsize=(3.8,3.2)) #4.5,4
   sns.lineplot(data=frap_results[frap_results['BG'] == 0],x='Time',y='scaled_intensity',hue='condition', 
             marker='o', ci='sd',linestyle='', err_style="bars", 
             hue_order = conditions,palette=colours) 
   plt.legend(title='',frameon=False,labels = 
           [x+' (N = ' + str(len(set(frap_results[frap_results['condition'] == x]['FRAP'])))+')' for x in conditions])
   plt.ylabel('Scaled intensity',fontsize=14)
   plt.xlabel('Time (s)',fontsize=14)
   plt.xlim(xmin=-2,xmax=maxT)
   plt.ylim([0,1.1])
   plt.tight_layout()
   plt.savefig('{}/{}.pdf'.format(outdir,figname),transparent=True)

def plot_curves_son(frap_results,conditions,colours,maxT,figname,outdir):
   ax_scale = 3.333 
   f, (ax1, ax2) = plt.subplots(ncols=1, nrows=2,
                             sharex=True,figsize=(3.8,3.5),gridspec_kw={'height_ratios': [1, ax_scale]})
   ax1 = sns.lineplot(data=frap_results[frap_results['BG'] == 0],x='Time',y='scaled_intensity',hue='condition',
             marker='o', ci='sd',linestyle='', err_style="bars",
             hue_order = conditions,palette=colours,ax=ax1)
   ax2 = sns.lineplot(data=frap_results[frap_results['BG'] == 0],x='Time',y='scaled_intensity',hue='condition',
             marker='o', ci='sd',linestyle='', err_style="bars",
             hue_order = conditions,palette=colours,ax=ax2)

   ax1.set_ylim(0.9, 1.05)
   ax2.set_ylim(0, 0.5)
   ax1.get_xaxis().set_visible(False)
   ax1.set_ylabel("")
   ax2.set_ylabel("")
   f.text(-0.03, 0.45, "Scaled intensity", fontsize = 14, va="center", rotation="vertical")


   ax1.spines['bottom'].set_visible(False)
   ax2.spines['top'].set_visible(False)
   d = .01  # how big to make the diagonal lines in axes coordinates
   kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False)
   ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
   ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
   kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
   ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
   ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

   ax2.set_xlim(-2,maxT)
   ax2.set_xlabel('Time (s)',fontsize=14)

   ax1.get_legend().remove()
   ax2.get_legend().remove()
   ax2.legend(loc=(1.025, 0.5), title="",frameon=False, labels =
           [x+' (N = ' + str(len(set(frap_results[frap_results['condition'] == x]['FRAP'])))+')' for x in conditions])
   plt.subplots_adjust(hspace=0.1)
   #plt.tight_layout()
   plt.savefig('{}/{}.pdf'.format(outdir,figname),transparent=True,bbox_inches='tight')

def import_data(frap_results,conditions):
   #print(conditions)
   frap_results['condensate_by_condition'] = frap_results['condition'] + ' ' + pd.Series([str(int(i)) for i in frap_results['FRAP']])
   frap_results['scaled_intensity'] = 0

   cond, recovery = [],[]
   for condensate in frap_results['condensate_by_condition'].dropna().unique().tolist():
      subdf = frap_results[frap_results['condensate_by_condition'] == condensate]
      cd = subdf['condition'].tolist()[0]
      if cd in conditions:
         cond.append(cd)
         max_in = np.max(subdf['Mean'])
         min_in = np.mean(subdf[subdf['BG'] == 1]['Mean'])
         subsub = subdf[subdf['BG'] == 0]
         scaled_i = (subdf['Mean']-min_in)/(max_in-min_in)
         frap_results.loc[frap_results['condensate_by_condition'] == condensate,'scaled_intensity'] = scaled_i    
         recovery.append(np.nanmedian(scaled_i[-10:]))
   recovery_df = pd.DataFrame(list(zip(cond,recovery)),
               columns =['condition', 'recovery'])
   for condition in conditions:
      print('mean recovery',condition,recovery_df[recovery_df['condition'] == condition]['recovery'].mean())
   return(frap_results,recovery_df)

def main():
   from argparse import ArgumentParser
   parser = ArgumentParser(description='Plot violin plot of recoveries for individual sites + recovery curves summarized per condition for FRAP')
   parser.add_argument('-x','--xlsx',type=str,required=True,help='xlsx file with FRAP results (must have columns labelled "Mean","Time","condition","BG", and "FRAP" at a minimum')
   parser.add_argument('-s','--conditions',nargs='+',type=str,required=False,help='list of subset of conditions from xlsx to plot')
   parser.add_argument('-c','--colours',nargs='+',type=str,required=False,help='list of colours for plot')
   parser.add_argument('-m','--maxT',type=int,required=False,help='max time to include for recovery curve')
   parser.add_argument('-n','--name',type=str,required=False,help='name for output plots')

   args = parser.parse_args()
   timedf = pd.read_excel(args.xlsx)
   if args.conditions is None:
       print(timedf['condition'])
       conditions = sorted(list(set(timedf['condition'].dropna().tolist())))
   else:
       conditions = args.conditions
   if not args.colours:
       colours = sns.color_palette("Paired", len(conditions))
   else:
       colours = args.colours
   if not args.maxT:
       maxT = max(timedf['Time']+1)
   else:
       maxT = args.maxT
   if not args.name:
       name = "frap"
   else:
       name = args.name
   
   outd = '' #set output directory
   timedf,recoveries = import_data(timedf,conditions)
   plot_curves(timedf,conditions,colours,maxT,name,outd) #standard
   #plot_curves_son(timedf,conditions,colours,maxT,name,outd) #to plot low-recovery curves for SON, with broken y-axis

if __name__ == "__main__":
   main()
