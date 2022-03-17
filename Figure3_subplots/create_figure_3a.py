import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

'''Gather all junctions/types in pickled output across all strains'''
junction_type = "Normal" #"Inverted"

currDir = os.path.dirname(os.path.realpath(__file__))
parent_path = os.path.abspath(os.path.join(currDir, os.pardir))
folder_name = '/pikl_pipeline_output_petites/'
path = parent_path + folder_name
#ori1, ori2, ori3, ...,ori8
origin_locations_forlines = [[4312, 4012], [32231, 32501], [54840, 54567], [56832, 56567], [82329, 82600] , [45227, 47927], [30220, 30594], [12780, 12510]]

allsamples_normal_junctions = []
allsamples_inverted_junctions = []

for file in os.listdir(path):
    with open(path + file, 'rb') as f:
        per_sample_pikl = pickle.load(f)

        primary_aln_key_count = per_sample_pikl[5]  # [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]
        alternate_aln_key_count = per_sample_pikl[6]  # [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]
        junction_ref_boundaries_dict = per_sample_pikl[
            7]  # key: abs(junction num), value: [[LL, LH, L mean], [HL, HH, H mean], [LL std, LH std, L mean std], [HL std, HH std, H mean std]
        junction_type_dict = per_sample_pikl[
            8]  # key: abs(junction num), value: True or False, (True meaning inverted junction type, False, noninverted)
        junction_counts_dict = per_sample_pikl[
            9]  # key : abs(junction num), value: count (this is for all true junctions whether or not they show up in repeats)

        #populate all samples lists
        for key in junction_ref_boundaries_dict.keys():
            if junction_type_dict[key]=='False':
                if key in primary_aln_key_count[0]:
                    allsamples_normal_junctions.append(
                        [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]])

                else:
                    for ele in alternate_aln_key_count:
                        if key in ele[0]:
                            allsamples_normal_junctions.append(
                                [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]])

            else:
                if key in primary_aln_key_count[0]:
                    allsamples_inverted_junctions.append(
                        [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]])

                else:
                    for ele in alternate_aln_key_count:
                        if key in ele[0]:
                            allsamples_inverted_junctions.append(
                                [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]])


# allsamples_normal_junctions.remove([1, 85777])
# allsamples_normal_junctions.remove([1, 85778])
# allsamples_normal_junctions.remove([1, 85776])
# allsamples_normal_junctions.remove([20, 85770])
# allsamples_normal_junctions.remove([48, 85750])
'''Plotting '''

size = 30
#plt.title('chrM breakpoint junctions \n (DBSCAN estimated number of clusters: %d)' % (n_clusters_normal + n_clusters_inverted), fontsize=size)

#plt.legend(ncol=1, loc='upper left', fontsize=14,handles=[fill_inv, fill_normal, class_member_inverted,class_member_normal, noise_member_inverted, noise_member_normal], bbox_to_anchor=(1.2,0.5), labels=["Inverted Junctions", "Normal Junctions", "Inverted junction class member","Normal junction class member", "Purported inverted junction noise", "Purported normal junction noise"])
toplot_normal = np.array(allsamples_normal_junctions)
toplot_inverted = np.array(allsamples_inverted_junctions)

if junction_type=="Normal":
    toplot = toplot_normal  #set what in particular to plot
else:
    toplot = toplot_inverted

marker='.'

plt.plot(toplot[:, 0], toplot[:, 1], marker, markerfacecolor='k',
                     markeredgecolor= 'k', markersize=15, alpha=0.5,  markeredgewidth=2)

plt.plot(toplot[:, 1], toplot[:, 0], marker, markerfacecolor='k',
                     markeredgecolor= 'k', markersize=15, alpha=0.5,  markeredgewidth=2)

line = np.arange(85779)

#plot all horizontal blocks corresponding to origins of replication
for ori in origin_locations_forlines:
    if ori[0] < ori[1]:
        fill_normal = plt.fill_between(line, [ori[0]]*len(line),[ori[1]]*len(line) , color='blue', alpha= 0.3)
    else:
        fill_normal = plt.fill_between(line, [ori[0]] * len(line), [ori[1]] * len(line), color='green', alpha=0.3)

#plot all vertical blocks corresponding to origins of replication
for ori in origin_locations_forlines:
    if ori[0] < ori[1]:
        fill_normal = plt.fill_between(np.arange(ori[0], ori[1]), 0, 85779 , color='blue', alpha= 0.3)
    else:
        fill_normal = plt.fill_between(np.arange(ori[1], ori[0]), 0, 85779, color='green', alpha=0.3)

#fill_inv = plt.fill_between(line, line, [max(line)]*len(line), color='green', alpha= 0.2)
plt.xlim(-4000, 85779)
plt.ylim(-4000, 85779)
plt.xlabel("chrM reference position (kbp)", fontsize=size)
plt.ylabel("chrM reference position (kbp)", fontsize=size)

# if junction_type== "Inverted":
#
#     plt.title("Inverted junction and replication origin locations", fontsize=size)
#
# else:
#
#     plt.title("Non-inverted junction and replication origin locations", fontsize=size)

plt.xticks([0, 20000, 40000, 60000, 80000], [0, 20, 40, 60, 80], fontsize=size)
plt.yticks([0, 20000, 40000, 60000, 80000], [0, 20, 40, 60, 80], fontsize=size)

'''Plotting origin locations below axes'''

#ax.margins(0.1)

line_pos = -1000

text_pos = -3800
text_pos_lower = -3000


#ax1.set_title("re")

#ax1.add_line(Line2D((xmin, xmax), (ymin*0.5, ymin*0.5), color='black', linewidth=4))
#plotting features of potential interest in its own subplot
feature_linewidth = 10
my_spline_width = 2

feature_size = 15

#plot false splines
plt.plot((0, 85779), (0,0), 'k-', linewidth = my_spline_width)
plt.plot((0, 0), (0,85779), 'k-', linewidth = my_spline_width)

#plot x axis origin locations
plt.plot((82329, 82600), (line_pos,line_pos), 'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(82496, text_pos, "ORI5", size = feature_size, ha='center')

plt.plot((56567, 56832), (line_pos,line_pos), 'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(56699, text_pos, "ORI4", size = feature_size, ha='left')

plt.plot((54567, 54840), (line_pos,line_pos), 'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(54703, text_pos, "ORI3", size = feature_size, ha='right')

plt.plot((45227, 47927), (line_pos,line_pos), 'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(46577, text_pos, "ORI6", size = feature_size, ha='center')

plt.plot((32231, 32501), (line_pos,line_pos), 'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(32411, text_pos, "ORI2", size = feature_size, ha='left')

plt.plot((30220, 30594), (line_pos,line_pos), 'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(30407, text_pos, "ORI7", size = feature_size, ha='right')

plt.plot((12510, 12780), (line_pos,line_pos), 'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(12645, text_pos, "ORI8", size = feature_size, ha='center')

plt.plot((4012, 4312), (line_pos,line_pos), 'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(4162, text_pos, "ORI1", size = feature_size, ha='center')

#plot y axis origin locations
plt.plot((line_pos,line_pos), (82329, 82600),  'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos, 82496, "ORI5", size = feature_size, ha='left', rotation=90)

plt.plot((line_pos,line_pos), (56567, 56832),  'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos, 56699,  "ORI4", size = feature_size, ha='left', va='bottom', rotation=90)

plt.plot((line_pos,line_pos),(54567, 54840),  'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos,54703,  "ORI3", size = feature_size, ha='left', va='top', rotation=90)

plt.plot((line_pos,line_pos),(45227, 47927),  'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos,46577,  "ORI6", size = feature_size, ha='left', rotation=90)

plt.plot((line_pos,line_pos),(32231, 32501),  'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos,32411,  "ORI2", size = feature_size, ha='left', va='bottom', rotation=90)

plt.plot((line_pos,line_pos),(30220, 30594),  'b-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos,30407,  "ORI7", size = feature_size, ha='left', va='top', rotation=90)

plt.plot((line_pos,line_pos),(12510, 12780),  'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos,12645,  "ORI8", size = feature_size, ha='left', rotation=90)

plt.plot((line_pos,line_pos),(4012, 4312),  'g-', linewidth = feature_linewidth, solid_capstyle="butt")
plt.text(text_pos, 4162,  "ORI1", size = feature_size, ha='left', rotation=90)

#legend creation
green_patch = mpatches.Patch(color='green', label='Negative strand origins', alpha=0.3)
blue_patch = mpatches.Patch(color='blue', label='Positive strand origins', alpha=0.3)
#plt.legend(ncol=1, handles=[green_patch, blue_patch], loc='upper left', bbox_to_anchor=(1.0,0.5), fontsize=size-5)
plt.legend(ncol=1, handles=[green_patch, blue_patch], loc='upper left', bbox_to_anchor=(0.5,0.3), fontsize=size-10)
#plt.legend(ncol=1, handles=[green_patch, blue_patch], loc='best', fontsize=size-5)

figure = plt.gcf()
figure.set_size_inches(12,10)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir +os.sep

folder_name = os.sep + "junction_dot_plot" + os.sep
path = currDir + folder_name

print path
if not os.path.exists(path):
    os.makedirs(path)
    print path

#plt.savefig(path +  "normal_junctions", bbox_inches="tight", DPI=300)
if junction_type == "Inverted":
    plt.savefig(path +  "inverted_junctions", bbox_inches="tight", dpi=600)
else:
    plt.savefig(path + "normal_junctions", bbox_inches="tight", dpi=600)









