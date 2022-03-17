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

petite_frac = {}
petite_frac_confident = {}
petite_frac_contained = {}

petite_content_confident = {}
WT_content_confident = {}

with open(folder + "petites_Petite_like_fraction.txt", 'r') as petite:
	for line in petite:
		name = line.strip().split()[0].split('_')[0]
		petite_frac[name] = float(line.strip().split()[1])
		petite_frac_confident[name]=float(line.strip().split()[2])
		petite_frac_contained[name]=float(line.strip().split()[3])

		petite_content_confident[name] = int(line.strip().split()[4])
		WT_content_confident[name] = int(line.strip().split()[5])

# fam1 = ["2I1", "21I1", "22I1"] #family 1a
# fam2 = ["18I1", "20I1", "10I1", "4I1", "5I1"] #family 2
# fam3 = ["3I1", "13I1", "19I1", "17I1", "23I1", "8I1", "6I1"] #family 3
# fam4 = ["9I1", "14I1", "24I1"] #family 4a
# fam5= ["1I1", "15I1", "16I1"] #family 4b
# fam6 = ["13I2", "20I2", "21I2"] #family 5
# fam7 = ["7I2", "8I2", "9I2"] #family 6
# fam8 = ["5I2", "22I2", "6I2"] #family 7
# fam9 = ["10I2", "11I2"] #family 8a
# fam10 = ["23I2", "24I2"] #family 9

# contingency_1 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam1]
# contingency_2 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam2]
# contingency_3 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam3]
# contingency_4 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam4]
# contingency_5 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam5]
# contingency_6 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam6]
# contingency_7 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam7]
# contingency_8 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam8]
# contingency_9 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam9]
# contingency_10 = [[WT_content_confident[fam_member],petite_content_confident[fam_member]] for fam_member in fam10]

####
# contingencies = [contingency_1, contingency_2, contingency_3, contingency_4, contingency_5, contingency_6, contingency_7, contingency_8, contingency_9, contingency_10]
# pvalue_allgroups = [chi2_contingency(np.array(contingency))[1] for contingency in contingencies]

#rearranging into high/low variance groups (arbitrary)
fam1 = ["18I1", "20I1", "10I1", "4I1", "5I1"] #family 2
fam2 = ["1I1", "15I1", "16I1"] #family 4b
fam3 = ["2I1", "21I1", "22I1"] #family 1a
fam4 = ["9I1", "14I1", "24I1"] #family 4a
fam5 = ["7I2", "8I2", "9I2"] #family 6
fam6 = ["5I2", "22I2", "6I2"] #family 7
fam7 = ["3I1", "13I1", "19I1", "17I1", "23I1", "8I1", "6I1"] #family 3
fam8 = ["23I2", "24I2"] #family 9
fam9 = ["10I2", "11I2"] #family 8a
fam10 = ["13I2", "20I2", "21I2"] #family 5

#family_colours = ["#2CEE0E", "#F3715A", "#E3D200", "#5C2D91", "#ED1C24", "#72BF44", "#00FFFF", "#9D85BE", "#0369A3"]

# palette = sns.color_palette(["#2CEE0E", "#F3715A", "#E3D200", "#5C2D91", "#5C2D91", "#ED1C24", "#72BF44", "#00FFFF", "#9D85BE", "#0369A3"])

#rearranging to high/low
palette = sns.color_palette(["#F3715A", "#5C2D91", "#2CEE0E", "#5C2D91", "#72BF44", "#00FFFF", "#E3D200", "#0369A3", "#9D85BE", "#ED1C24"])

family_sample_names = [fam1, fam2, fam3, fam4, fam5, fam6, fam7, fam8, fam9, fam10]

family_mt_fracs = [[petite_frac[key]for key in fam1], [petite_frac[key]for key in fam2], [petite_frac[key]for key in fam3], [petite_frac[key]for key in fam4],
					[petite_frac[key]for key in fam5], [petite_frac[key]for key in fam6], [petite_frac[key]for key in fam7], [petite_frac[key]for key in fam8],
					[petite_frac[key]for key in fam9], [petite_frac[key]for key in fam10]]

family_mt_fracs_confident = [[petite_frac_confident[key]for key in fam1], [petite_frac_confident[key]for key in fam2], [petite_frac_confident[key]for key in fam3], [petite_frac_confident[key]for key in fam4],
					[petite_frac_confident[key]for key in fam5], [petite_frac_confident[key]for key in fam6], [petite_frac_confident[key]for key in fam7], [petite_frac_confident[key]for key in fam8],
					[petite_frac_confident[key]for key in fam9], [petite_frac_confident[key]for key in fam10]]

family_mt_fracs_contained = [[petite_frac_contained[key]for key in fam1], [petite_frac_contained[key]for key in fam2], [petite_frac_contained[key]for key in fam3], [petite_frac_contained[key]for key in fam4],
					[petite_frac_contained[key]for key in fam5], [petite_frac_contained[key]for key in fam6], [petite_frac_contained[key]for key in fam7], [petite_frac_contained[key]for key in fam8],
					[petite_frac_contained[key]for key in fam9], [petite_frac_contained[key]for key in fam10]]

print
petite_version = family_mt_fracs_confident

#now do jitter plot
size=32
means = [np.mean(ele) for ele in petite_version]
stds = [np.std(ele) for ele in petite_version]

sorted_stds, sorted_order = (list(t) for t in zip(*sorted(zip(stds, [1,2,3,4,5,6,7,8,9,10]), reverse=False)))

print sorted_stds
print sorted_order
#do sorting for plotting based on mean frequency (to look pretty only)
#test_means, sorted_family_sample_names = (list(t) for t in zip(*sorted(zip(means, family_sample_names), reverse=True)))
test_means, sorted_family_sample_names = means, family_sample_names

print sorted_family_sample_names

sorted_family_mt_fracs_mc = []
sorted_means = []
names_recorded = []

for mean in test_means:
	for names, ele  in zip(family_sample_names, petite_version):
		if mean == np.mean(ele) and names not in names_recorded:
			#print ele, mean
			sorted_family_mt_fracs_mc.append(ele)
			sorted_means.append(mean)
			names_recorded.append(names)

#print sorted_family_mt_fracs_mc
labels = []
output_fam = []
# test_labels = [r'$1_a$', '2', '3', r'$4_a$', r'$4_b$', '5', '6', '7', r'$8_a$', '9']

test_labels = [ '2', r'$4_b$',r'$1_a$', r'$4_a$','6', '7', '3','9',  r'$8_a$', '5']

for i, fam in enumerate(sorted_family_mt_fracs_mc):
	labels.extend([test_labels[i]]*len(fam))
	output_fam.append(fam)

print output_fam
toplot = []
for fam in sorted_family_mt_fracs_mc:
	toplot.extend(fam)

d = {'Sample': labels, 'Alternate mt content fraction' : toplot}
df = pd.DataFrame(d)

#sns.violinplot(data = df, x='Junction', y='Distance to closest origin (bp)', color="0.8", cut=0, inner=None)
sns.stripplot(data = df, x='Sample', y='Alternate mt content fraction', jitter=True, zorder=1, size=10, palette=palette)
#plt.scatter(x=range(len(sorted_means)), y=sorted_means, c='k')

#data = [sorted(np.random.normal(0, std, 10)) for std in range(1, 5)]
#plt.title("Alternate structure content fractions", fontsize=size)

plt.ylabel("Alternate structure\nbp fraction", fontsize=size)
#plt.ylim(-0.00001, 0.3)
#plt.yscale('symlog')
plt.xlabel("Petite family/sub-family", fontsize=size)
plt.xticks(fontsize=size-5, rotation=0)
plt.yticks(fontsize=size, rotation=0)
#plt.ylim(-0.005, 0.02)

# plt.text(0, 0.10, r"$\chi^2$ p-value=" + "\n%.2E" %(pvalue_allgroups[0]), fontsize=size, ha='center')
# plt.text(1, 0.10, r"$\chi^2$ p-value=" + "\n%.2E" %(pvalue_allgroups[1]), fontsize=size, ha='center')
# plt.text(2, 0.05, r"$\chi^2$ p-value=" + "\n%.2E" %(pvalue_allgroups[2]), fontsize=size, ha='center')
# plt.text(3, 0.05, r"$\chi^2$ p-value=" + "\n%.2E" %(pvalue_allgroups[3]), fontsize=size, ha='center')
# plt.text(4, 0.10, r"$\chi^2$ p-value=" + "\n%.2E" %(pvalue_allgroups[4]), fontsize=size, ha='center')
# plt.text(5, 0.12, r"$\chi^2$ p-value=" + "\n%.2E" %(pvalue_allgroups[5]), fontsize=size, ha='center')

#plt.yscale('log')

figure = plt.gcf()
figure.set_size_inches(12, 10)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

folder_name = os.sep + "nature_of_diversity" + os.sep
path = currDir + folder_name

if not os.path.exists(path):
    os.makedirs(path)

#save histogram figure
plt.savefig(path + "fraction_alternate_content_all_families_nostructure.png", bbox_inches="tight", dpi=300, format='PNG')

