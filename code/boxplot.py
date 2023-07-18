##### For creating geneomic bins distribution in the form of Boxplots ######

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from statannot import add_stat_annotation
import seaborn as sns; sns.set_theme(color_codes=True)
import seaborn as sns; sns.set_theme(color_codes=True)
from seaborn import lmplot
from sklearn import preprocessing
import sys
import seaborn as sns
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

df = pd.read_csv('H3K27ac_input.csv')


dd=pd.melt(df,id_vars=['Group'],value_vars=['HCT116','RH4'],var_name='Gene Regulation')
ax = sns.boxplot(x='Group',y='value',data=dd,hue='Gene Regulation')


sns.set_style('whitegrid')
ax.set(ylim=(-0.5, 1.5))
ax.set_ylabel("Zscore Normalized Read Count",fontsize=25, fontweight="bold")
ax.set_xlabel("Genomic Region",fontsize=25, fontweight="bold")

ax.figure.set_size_inches(25,8)
#ax.get_legend().remove()            ### If don't want labels on the image
#ax.legend(loc=1)                   ### '1' will shift label on top, default is down.

plt.savefig('Output_image', dpi=600)
