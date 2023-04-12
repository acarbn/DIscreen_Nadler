import matplotlib.pyplot as plt
import seaborn as sns 
import matplotlib.cm as cm
import pandas as pd
import numpy as np
from scipy import stats

# set styles
sns.set_style("ticks", rc={'font.family':'sans-serif', 'font.sans-serif':'Droid Sans'})
plt.style.use('./utils/domain_ins.mplstyle')
plt.rcParams['svg.fonttype'] = 'none'


def r(x, y):
    return stats.pearsonr(x, y)[0]

def create_enrichment_fig(data_norm_1, data_norm_2, combination, condition, prot_dict, out_folder):
    #plt.clf()
    #print(data_norm_1, data_norm_2)
    plt.figure(figsize=(10,3))
    sns.set_context(rc = {'patch.linewidth': 0.0}) 
    try:
        clrs_1 = ['grey' if x == -10 else '#87001D' if pd.isnull(x) else '#00B8B8' for x in data_norm_1.log]
        clrs_2 = ['grey' if x == -10 else '#87001D' if pd.isnull(x) else '#174950' for x in data_norm_2.log]
        clrs_1[0] = 'white'
        clrs_2[0] = 'white'
        ax = sns.barplot(data=data_norm_1.fillna(-10), y='log', x='position', palette=clrs_1, label='Rep-1', alpha=.8)
        ax2 = sns.barplot(data=data_norm_2.fillna(-10), y='log', x='position', palette=clrs_2, label='Rep-2', alpha=.8)
        labels=data_norm_1.loc[::20,'position']
        ax.legend([])
        ax2.legend([])
        for patch in ax2.patches:
            patch.set_width(1)
    except:
        try:
            clrs_1 = ['grey' if x == -10 else '#87001D' if pd.isnull(x) else '#00B8B8' for x in data_norm_1.log]
            clrs_1[0] = 'white'
            ax = sns.barplot(data=data_norm_1.fillna(-10), y='log', x='position', palette=clrs_1, label='Rep-1')
            labels=data_norm_1.loc[::20, 'position']
        except:
            clrs_2 = ['grey' if x == -10 else '#87001D' if pd.isnull(x) else '#174950' for x in data_norm_2.log]
            clrs_2[0] = 'white'
            ax = sns.barplot(data=data_norm_2.fillna(-10), y='log', x='position', palette=clrs_2, label='Rep-2')
            labels=data_norm_2.loc[::20, 'position']
    plt.title(f"{combination.replace('_', ' ')} Enrichment: {condition[0][1:]}")
    ax.set_ylabel("Log2 enriched read counts")
    plt.xlabel("Insertion site")
    for patch in ax.patches:
        patch.set_width(1)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.yaxis.set_tick_params(width=2, which='both')
    ax.xaxis.set_tick_params(width=2)
    plt.xticks(ticks= np.arange(0, len(prot_dict[combination.split('_')[0]]), 20.0), labels=labels)
    plt.xlim(left=-0.5, right=len(prot_dict[combination.split('_')[0]]))
    plt.ylim(bottom=-10)
    sns.despine(top=True, right=True)
    plt.savefig(f"{out_folder}/figures/enrichment_{combination}_{condition[0][1:]}.svg")
    plt.show()
    plt.close()


def create_fold_change_fig(data, misc, combination, condition, prot_dict, out_folder):
    #plt.clf()
    plt.figure(figsize=(10,3.5))
    sns.set_context(rc = {'patch.linewidth': 0.0})
    clrs = ['grey' if (x+1 in misc) else '#E60234' for x in range(len(data))]
    ax = sns.barplot(data.index, data.values, palette=clrs, label='Rep-1')
    labels=data.index[::20]
    plt.title(f"Switching fold changes: {combination.replace('_', ' ')} Enrichment: {condition[0][1:].split('_')[0]}")
    ax.set_ylabel("Fold change induced/uninduced")
    plt.xlabel("Insertion site")
    for patch in ax.patches:
        patch.set_width(1)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.yaxis.set_tick_params(width=2, which='both')
    ax.xaxis.set_tick_params(width=2)
    plt.xticks(ticks= np.arange(0, len(prot_dict[combination.split('_')[0]]), 20.0), labels=labels)
    plt.xlim(left=-0.5, right=len(prot_dict[combination.split('_')[0]]))
    plt.ylim(bottom=-50, top=50)
    sns.despine(top=True, right=True)
    plt.savefig(f"{out_folder}/figures/fold_change_{combination}_{condition[0][1:]}.svg")

def correlation_plot(data, combination, condition, out_folder):
    plt.figure(figsize=(5,5))
    plt.rcParams['axes.linewidth'] = 2
    ax = sns.regplot(data=data, x='Replicate-2', y='Replicate-1', color='grey', 
        scatter_kws={"color": 'grey', 'alpha':.4, 'linewidth':0}, line_kws={"color": 'grey', 
        'label':f"Pearson's r: {round(r(data['Replicate-2'], data['Replicate-1']), 2)}"})
    plt.xlabel("Rep-2 log2 enriched read counts")
    plt.ylabel("Rep-1 log2 enriched read counts")
    plt.xlim(data.min().min()-1, data.max().max() +1)
    plt.ylim(data.min().min()-1, data.max().max() +1)
    sns.despine()
    plt.legend()
    plt.title(f"{combination.replace('_', ' ')} Enrichment: {condition[0][1:]}")
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    plt.savefig(f"{out_folder}/figures/correlation_{combination}_{condition[0][1:]}.svg")
    plt.show()
    plt.close()

def color_map_color(value, cmap_name='viridis', vmin=0, vmax=1):
    '''
    color map for protein depiction
    '''
    # norm = plt.Normalize(vmin, vmax)
    norm = plt.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(cmap_name)  # PiYG
    rgb = cmap(norm(value))[:3]  # will return rgba, we take only first 3 so we get rgb
    color = plt.colors.rgb2hex(rgb)
    return color

'''
def correlation_plot(data, combination, property):
    data = data.loc[data[property].notna()]
    plt.figure(figsize=(5,5))
    plt.rcParams['axes.linewidth'] = 2
    g = sns.lmplot(data=data, x=property, y='enrichment', hue='variable', palette={'norm':'#008080'}, ci=None, 
        scatter_kws={'alpha':.3, 'linewidth':0}, line_kws={'alpha':1})
    plt.xlabel(property)
    plt.ylabel("Log2 variant enrichment")
    sns.despine()
    g._legend.remove()
    plt.title(f"Correlation between enrichment and {property}: {combination}")
    for ax in g.axes.flat:
        ax.yaxis.set_tick_params(width=2)
        ax.xaxis.set_tick_params(width=2)
    #plt.savefig(f"{in_folder}/coverage_plots/log_transformed/correlation_{combination}_{condition[0][1:]}.png")
'''
def swarm_plot(data, combination, property):
    data = data.loc[data[property].notna()]
    plt.figure(figsize=(5,5))
    plt.rcParams['axes.linewidth'] = 2
    ax = sns.boxplot(data=data, x=property, y='enrichment', linewidth=2, color='white', fliersize=0)
    ax = sns.swarmplot(data=data, x=property, y='enrichment', alpha=.5, linewidth=0, color='grey')
    plt.setp(ax.lines, color="0")
    plt.xlabel(property)
    plt.ylabel("Log2 variant enrichment")
    sns.despine()
    plt.title(f"Correlation between enrichment and {property}: {combination}")
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    #plt.savefig(f"{in_folder}/coverage_plots/log_transformed/correlation_{combination}_{condition[0][1:]}.png")

map = sns.diverging_palette(170, 260, s=60, as_cmap=True)

def pairwise_correlation(data, name, out_folder):
    sns.set(font_scale = 2)
    
    data = data.corr(method='pearson')
    plt.figure(figsize=(10,10))
    plt.rcParams['axes.linewidth'] = 2
    ax = sns.heatmap(data=data, cmap='mako', cbar_kws={'label': "pearson's r"}, square=True)
    plt.title(f"Correlation between enrichments: {name}")
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    plt.savefig(f"{out_folder}/figures/{name}_correlation-1.svg")
'''
def pairwise_correlation(data, name, method = 'pearson'):
    sns.set_style("ticks")
    sns.set(font_scale = 1)
    data = data.corr(method=method)
    plt.figure(figsize=(10,10))
    ax = sns.heatmap(data=data, cmap=map, cbar_kws={'label': f"{method}'s r"}, square=True, vmin=-1, vmax=1)
    plt.title(f"Correlation between enrichments: {name}")
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    #plt.savefig(f"{in_folder}/coverage_plots/correlations/feature_correlation_{name}.png", bbox_inches="tight")
    plt.show()
    plt.close()
'''
def create_switch_fig(data, out_folder):
    data = data.fillna(0)
    plt.figure(figsize=(10,3)) 
    clrs_l = []
    clrs_d = []
    for idx, i  in data.iterrows():
        if i['light'] != 0:
            clrs_l.append('light')
        else:
            clrs_l.append('0')
        if i['dark'] != 0:
            clrs_d.append('dark')
        else:
            clrs_d.append('0')
    ax = sns.barplot(data=data, y=data['dark']-data['light'], x='position', color='black', edgecolor='black', bottom=data['light'], label='induced more active', alpha=1, zorder=0, linewidth=0)
    ax3 = sns.scatterplot(data=data, y='light', x=data['position']-1, hue=clrs_l, palette={'light':'#2582FF', '0':'#FF00'}, legend=False, edgecolor='none', s=50, zorder=2)
    ax4 = sns.scatterplot(data=data, y='dark', x=data['position']-1, hue=clrs_d, palette={'dark':'#808080', '0':'#FF00'}, legend=False, edgecolor='none', s=50, zorder=3)
    labels=data.position[::20]
    for patch in ax.patches:
        patch.set_width(2)
        patch.set_x(patch.get_x() - .5)
    
    plt.title('Light-switchable AraC-LOV hybrids')
    ax.set_ylabel("Log2 enriched read counts")
    plt.xlabel("Insertion site")
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.yaxis.set_tick_params(width=2, which='both')
    ax.xaxis.set_tick_params(width=2)
    plt.xticks(ticks= np.arange(0, 292, 20.0), labels=labels)
    #plt.xlim(left=100, right=120)
    plt.ylim(bottom=-7.5, top=4)
    sns.despine(top=True, right=True)
    plt.savefig(f"{out_folder}/figures/AraC_switch-1.svg")