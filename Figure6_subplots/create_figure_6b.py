import matplotlib.pyplot as plt
import glob
import numpy as np
import pickle
import math
from scipy.stats import chi2_contingency
import seaborn as sns
import pandas as pd
import os
def rotate(l, n):
    return l[-n:] + l[:-n]

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

folder = path + os.sep + 'petite_fraction' + os.sep

YPD_frac = []
YPG_frac = []
YPD_frac_confident = []
YPG_frac_confident = []

petite_frac = []
petite_frac_confident = []
petite_frac_contained = []

with open(folder + "YPD_Petite_like_fraction.txt", 'r') as YPD:
	for line in YPD:
		YPD_frac.append(float(line.strip().split()[1]))
		YPD_frac_confident.append(float(line.strip().split()[2]))


with open(folder + "YPG_Petite_like_fraction.txt", 'r') as YPG:
	for line in YPG:
		YPG_frac.append(float(line.strip().split()[1]))
		YPG_frac_confident.append(float(line.strip().split()[2]))

with open(folder + "petites_Petite_like_fraction.txt", 'r') as petite:
	for line in petite:
		petite_frac.append(float(line.strip().split()[1]))
		petite_frac_confident.append(float(line.strip().split()[2]))
		petite_frac_contained.append(float(line.strip().split()[3]))

YPD_version = YPD_frac_confident
YPG_version = YPG_frac_confident
petite_version = petite_frac_confident

#now do jitter plot
size=32
means = [np.mean(YPD_version), np.mean(YPG_version), np.mean(petite_version)]

YPG_label = ['Grande YPG' for i in range(len(YPG_version))]
YPD_label = ['Grande YPD' for i in range(len(YPD_version))]
petite_label = ['Petite YPD' for i in range(len(petite_version))]

labels = []
labels.extend(YPD_label)
labels.extend(YPG_label)
labels.extend(petite_label)

toplot = []
toplot.extend(YPD_version)
toplot.extend(YPG_version)
toplot.extend(petite_version)

d = {'Sample': labels, 'Alternate mt content fraction' : toplot}
df = pd.DataFrame(d)

#sns.violinplot(data = df, x='Junction', y='Distance to closest origin (bp)', color="0.8", cut=0, inner=None)
sns.stripplot(data = df, x='Sample', y='Alternate mt content fraction', color=".25", jitter=True, zorder=1, size=10)
ax = sns.boxplot(data=df, x='Sample', y='Alternate mt content fraction',fliersize=0)
#plt.scatter(x=range(len(means)), y=means, c='k', s=50)

#data = [sorted(np.random.normal(0, std, 10)) for std in range(1, 5)]
#plt.title("Petite-like structure content fractions", fontsize=size)

plt.ylabel("Alternate structure\nbp fraction", fontsize=size)
plt.xlabel("Colony type/media conditions", fontsize=size)
plt.xticks(fontsize=size-5, rotation=0)
plt.yticks(fontsize=size, rotation=0)
#plt.ylim(-0.005, 0.02)

figure = plt.gcf()
figure.set_size_inches(12, 10)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

folder_name = os.sep + "nature_of_diversity" + os.sep
path = currDir + folder_name

if not os.path.exists(path):
    os.makedirs(path)

#save histogram figure
plt.savefig(path + "fraction_petite-like_content_confident.png", bbox_inches="tight", dpi=300, format='PNG')

