import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from matplotlib.pyplot import cm
import matplotlib.lines as mlines

currDir = os.path.dirname(os.path.realpath(__file__))
parent_path = os.path.abspath(os.path.join(currDir, os.pardir))

'''Must change for nuclear vs mt for annotated features filename and header, chr length etc'''
header = 6 # how many lines before data actually starts?
featurefile = "All Annotated Sequence Features-chrmt-1..85779_culled.sqn"
features = {}
raw_lines = []
line_feature_pos = []

def sequin_format(argument):
    if len(argument) != 0 and len(argument[0].strip().split()) <= 2:
        return argument[0].strip().split()[1]
    else:
        return "N/A"

with open(os.path.join(parent_path, featurefile), 'r') as annotated:
    # loop to get feature line positions in file
    for pos, line in enumerate(annotated):
        # if we are starting a new feature
        if not line.startswith("\t") and (pos >= header):
            line_feature_pos.append(pos)
        if (pos >= header):  # if we are outside the whole "chromosome" annotation which we definitely don't want
            raw_lines.append(line)

for i in range(0, len(line_feature_pos) - 1):
    startline = line_feature_pos[i] - header
    endline = line_feature_pos[i + 1] - header

    current = raw_lines[startline:endline]

    ID = [s for s in current if "id" in s]
    NAME = [s for s in current if "name" in s]
    NOTE = [s for s in current if "note" in s]
    GENE = [s for s in current[1:] if "gene" in s]
    START = current[0].strip().split()[0]
    STOP = current[0].strip().split()[1]
    GENETYPE = current[0].strip().split()[2]

    features[sequin_format(NAME)] = [GENETYPE, sequin_format(GENE), sequin_format(ID), sequin_format(NAME),
                                     NOTE[0].strip().split('\t')[-1], START, STOP]

# do the same for the last one:
current = raw_lines[endline:]

ID = [s for s in current if "id" in s]
NAME = [s for s in current if "name" in s]
NOTE = [s for s in current if "note" in s]
GENE = [s for s in current[1:] if "gene" in s]
START = current[0].strip().split()[0]
STOP = current[0].strip().split()[1]
GENETYPE = current[0].strip().split()[2]

features[sequin_format(NAME)] = [GENETYPE, sequin_format(GENE), sequin_format(ID), sequin_format(NAME),
                                 NOTE[0].strip().split('\t')[-1], START, STOP]

feature_bars = {}
for feature_key in features.keys():
    # replace gene with name for key in features_in_read dict if possible
    feature_gene = features[feature_key][1]
    feature_name = features[feature_key][3]
    feature_start = int(features[feature_key][5])
    feature_end = int(features[feature_key][6])

    if feature_gene != 'N/A':
        identifier = feature_gene
    else:
        identifier = feature_name

    feature_bars[identifier] = [min(feature_start, feature_end), max(feature_start, feature_end)]

print feature_bars

#construct colour map
feature_names = []
feature_locations = []
reference_starts = []

for key in feature_bars:
    feature_names.append(key)
    feature_locations.append(feature_bars[key])
    reference_starts.append(feature_bars[key][0])

temp = reference_starts
reference_starts, feature_names = (list(t) for t in zip(*sorted(zip(reference_starts, feature_names))))

temp, feature_locations = (list(t) for t in zip(*sorted(zip(temp, feature_locations))))

print feature_names
print feature_locations

n = len(feature_names)
colours = cm.prism(np.linspace(0,1,n))

'''Gather all junctions/types in pickled output across all strains'''
junction_type = "Normal" #"Inverted"

folder_name = '/pikl_pipeline_output_petites/'
path = parent_path+ folder_name
print(path)
#ori1, ori2, ori3, ...,ori8
families = ['_CP1', '_CP2', '_CP3', '_CP5', '_P2', '_P3', '_P4', '_P5', '_P7']
family_colours = ["#2CEE0E", "#F3715A", "#E3D200", "#5C2D91", "#ED1C24", "#72BF44", "#00FFFF", "#9D85BE", "#0369A3"]

origin_locations_forlines = [[4312, 4012], [32231, 32501], [54840, 54567], [56832, 56567], [82329, 82600] , [45227, 47927], [30220, 30594], [12780, 12510]]

allsamples_junctions_petites_normal = [] #all samples in Petites (in repeat structures)
allsamples_junctions_petites_inverted = []
allsamples_junctions_grandes_normal = [] #all samples in Grandes (in repeat structures if present, otherwise 1, 85779)
allsamples_colours_normal = []
allsamples_colours_inverted = []

junction_count_per_family = {}
for file in os.listdir(path):
    sample_name =  '_' + file.strip().split('/')[0].split('_')[1]
    test_family = ""
    print(file)
    with open(path + file, 'rb') as f:
        per_sample_pikl = pickle.load(f)

        primary_aln_key_count = per_sample_pikl[
            5]  # [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]
        alternate_aln_key_count = per_sample_pikl[
            6]  # [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]

        junction_ref_boundaries_dict = per_sample_pikl[
            7]  # key: abs(junction num), value: [[LL, LH, L mean], [HL, HH, H mean], [LL std, LH std, L mean std], [HL std, HH std, H mean std]
        junction_type_dict = per_sample_pikl[
            8]  # key: abs(junction num), value: True or False, (True meaning inverted junction type, False, noninverted)
        junction_counts_dict = per_sample_pikl[
            9]  # key : abs(junction num), value: count (this is for all true junctions whether or not they show up in repeats)

        keys_in_repeats_inverted = []
        keys_in_repeats_normal = []
        for key in junction_ref_boundaries_dict.keys():
            if junction_type_dict[key] == 'False':
                if key in primary_aln_key_count[0] and key not in keys_in_repeats_normal:
                    keys_in_repeats_normal.append(key)

                else:
                    for ele in alternate_aln_key_count:
                        if key in ele[0] and key not in keys_in_repeats_normal:
                            keys_in_repeats_normal.append(key)

            elif junction_type_dict[key] == 'True':
                if key in primary_aln_key_count[0] and key not in keys_in_repeats_inverted:
                    keys_in_repeats_inverted.append(key)

                else:
                    for ele in alternate_aln_key_count:
                        if key in ele[0] and key not in keys_in_repeats_inverted:
                            keys_in_repeats_inverted.append(key)

        #assign colours
        for i, family_name in enumerate(families):
            if family_name in sample_name:
                colour = family_colours[i]
                test_family =  family_name

        allsamples_colours_normal.extend([colour for key in keys_in_repeats_normal])
        allsamples_colours_inverted.extend([colour for key in keys_in_repeats_inverted])

        allsamples_junctions_petites_normal.extend([
            [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in
            keys_in_repeats_normal if abs(int(junction_ref_boundaries_dict[key][0][2]) - int(junction_ref_boundaries_dict[key][1][2])) <85700])

        test1 = [[junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in keys_in_repeats_normal if abs(int(junction_ref_boundaries_dict[key][0][2]) - int(junction_ref_boundaries_dict[key][1][2])) <85700]
        test2 = [[junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in keys_in_repeats_inverted]

        if test_family not in junction_count_per_family.keys():
            junction_count_per_family[test_family] = 0

        junction_count_per_family[test_family]+= len(test1)+len(test2)

        allsamples_junctions_petites_inverted.extend([
            [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in
            keys_in_repeats_inverted])


for key in junction_count_per_family.keys():
    print key, junction_count_per_family[key]

'''Plotting '''

size = 30
#plt.title('chrM breakpoint junctions \n (DBSCAN estimated number of clusters: %d)' % (n_clusters_normal + n_clusters_inverted), fontsize=size)

#plt.legend(ncol=1, loc='upper left', fontsize=14,handles=[fill_inv, fill_normal, class_member_inverted,class_member_normal, noise_member_inverted, noise_member_normal], bbox_to_anchor=(1.2,0.5), labels=["Inverted Junctions", "Normal Junctions", "Inverted junction class member","Normal junction class member", "Purported inverted junction noise", "Purported normal junction noise"])
petites_normal_index = len(allsamples_junctions_petites_normal) -1
toplot_petites = allsamples_junctions_petites_normal
toplot_petites.extend(allsamples_junctions_petites_inverted)

colours_for_plotting = allsamples_colours_normal
colours_for_plotting.extend(allsamples_colours_inverted)
#colours_for_plotting = np.array(colours_for_plotting)

print colours_for_plotting

toplot_grandes = [[0, 85779], [0, 85779], [0, 857779], [0,85779], [0, 85779], [0, 85779], [30325, 45316], [31234, 32996], [0, 85779], [0, 85779], [0, 85779], [0, 85779], [30643, 32589]]
toplot_grandes = [[30325, 45316], [31234, 32996], [30643, 32589]]

toplot_petites = np.array(toplot_petites)
toplot_grandes = np.array(toplot_grandes)

#toplot_inverted = np.array(allsamples_inverted_junctions)

fig, ax = plt.subplots()

print colours_for_plotting
petites_normal = ax.scatter(toplot_petites[:petites_normal_index+1, 0], toplot_petites[:petites_normal_index+1, 1], color=colours_for_plotting[:petites_normal_index+1], edgecolor='k', s=120, marker='o')
petites_inverted = ax.scatter(toplot_petites[petites_normal_index+1:, 0], toplot_petites[petites_normal_index+1:, 1], color=colours_for_plotting[petites_normal_index+1:], edgecolor='k', s=120, marker='^')
grandes_normal = ax.scatter(toplot_grandes[:, 1], toplot_grandes[:, 0], color=['black' for ele in toplot_grandes], edgecolor='k', s=120, marker='o')

print len(toplot_petites), "petite junction count"
print len(toplot_grandes), "grandes junction count"



plt.plot(np.arange(85779), np.arange(85779), 'k-', linewidth=2)
ax.set_facecolor('silver')

plt.text(5000, 80000, "Petite", fontsize=size)
plt.text(75000,5000,"Grande", fontsize=size)

line = np.arange(85779)

plt.xlim(-6000, 87779)
plt.ylim(-6000, 87779)
plt.xlabel("chrM reference position (kbp)", fontsize=size)
plt.ylabel("chrM reference position (kbp)", fontsize=size)

#plt.title("mtDNA junction locations from Nanopore sequencing of Petite/Grande samples", fontsize=size)

plt.xticks([0, 20000, 40000, 60000, 80000],[0,20,40,60,80], fontsize=size)
plt.yticks([0, 20000, 40000, 60000, 80000],[0,20,40,60,80], fontsize=size)

'''Plotting origin locations below axes'''

#ax.margins(0.1)

line_pos = -3000

text_pos_upper = -1800
text_pos_lower = -5000

#ax1.set_title("re")

#ax1.add_line(Line2D((xmin, xmax), (ymin*0.5, ymin*0.5), color='black', linewidth=4))
#plotting features of potential interest in its own subplot
feature_linewidth = 10
my_spline_width = 2

#plot false splines
plt.plot((0, 85779), (0,0), 'k-', linewidth = my_spline_width)
plt.plot((0, 0), (0,85779), 'k-', linewidth = my_spline_width)

#plot x axis values

feature_size = 20
upper=True
for i in range(len(feature_names)):
    if upper:
        plt.plot((feature_locations[i][0], feature_locations[i][1]), (line_pos, line_pos), color=colours[i], linewidth=feature_linewidth, solid_capstyle="butt")
        if feature_names[i] == 'ATP8' or feature_names[i] == 'ATP6' or feature_names[i] == 'ORI3':
            plt.text(np.mean((feature_locations[i][0], feature_locations[i][1])), text_pos_lower, feature_names[i],
                     size=feature_size, ha='right')
        else:
            plt.text(np.mean((feature_locations[i][0], feature_locations[i][1])), text_pos_lower, feature_names[i], size=feature_size, ha='center')
        upper=False

    else:

        plt.plot((feature_locations[i][0], feature_locations[i][1]), (line_pos, line_pos), color=colours[i],
                 linewidth=feature_linewidth, solid_capstyle="butt")
        if feature_names[i]=='ATP8' or feature_names[i]=='ATP6' or feature_names[i]=='ORI3':
            plt.text(np.mean((feature_locations[i][0], feature_locations[i][1])), text_pos_upper, feature_names[i],
                     size=feature_size, ha='right')
        else:
            plt.text(np.mean((feature_locations[i][0], feature_locations[i][1])), text_pos_upper, feature_names[i],
                     size=feature_size, ha='center')
        upper=True

#plot y axis origin locations
upper=False
for i in range(len(feature_names)):
    if upper:
        plt.plot((line_pos, line_pos),(feature_locations[i][0], feature_locations[i][1]), color=colours[i], linewidth=feature_linewidth, solid_capstyle="butt")
        if feature_names[i] == 'ATP8' or feature_names[i] == 'ATP6' or feature_names[i] == 'ORI3':
            plt.text(text_pos_lower-500, np.mean((feature_locations[i][0], feature_locations[i][1])), feature_names[i],
                     size=feature_size, ha='left', va='top', rotation=90)
        else:
            plt.text(text_pos_lower-500, np.mean((feature_locations[i][0], feature_locations[i][1])), feature_names[i],
                     size=feature_size, ha='left', va='center', rotation=90)
        upper=False

    else:
        plt.plot((line_pos, line_pos),(feature_locations[i][0], feature_locations[i][1]), color=colours[i],
                 linewidth=feature_linewidth, solid_capstyle="butt")
        if feature_names[i] == 'ATP8' or feature_names[i] == 'ATP6' or feature_names[i] == 'ORI3':
            plt.text(text_pos_upper, np.mean((feature_locations[i][0], feature_locations[i][1])), feature_names[i],
                     size=feature_size, ha='left', va='top', rotation=90)
        else:
            plt.text(text_pos_upper, np.mean((feature_locations[i][0], feature_locations[i][1])), feature_names[i],
                     size=feature_size, ha='left', va='center', rotation=90)
        upper=True


#add legend
class_member_inverted = mlines.Line2D([], [], marker='^', markerfacecolor='white',
             markeredgecolor='k', markersize=11, alpha=0.4, linestyle='None')

class_member_normal = mlines.Line2D([], [], marker='o', markerfacecolor='white',
             markeredgecolor='k', markersize=11, alpha=0.4, linestyle='None')

class_member_1 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[0],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_2 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[1],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_3 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[2],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_4 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[3],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_5 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[4],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_6 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[5],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_7 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[6],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_8 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[7],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_9 = mlines.Line2D([], [], marker='o', markerfacecolor=family_colours[8],
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')
class_member_grande = mlines.Line2D([], [], marker='o', markerfacecolor='black',
             markeredgecolor='k', markersize=11, alpha=1, linestyle='None')

#plt.legend(ncol=1, loc='center right', bbox_to_anchor=(1.0, 0.4), fontsize=size,handles=[class_member_inverted,class_member_normal], labels=["Inverted Junctions", "Non-inverted Junctions"])
plt.legend(ncol=1,  loc='center right', bbox_to_anchor=(1.005, 0.358), fontsize=size-6.5, handles=[class_member_inverted,class_member_normal, class_member_1, class_member_2, class_member_3, class_member_4, class_member_5, class_member_6, class_member_7, class_member_8, class_member_9, class_member_grande], labels=["Inverted breakpoints", "Non-inverted breakpoints","Petite family 1 (18 breakpoints/6 colonies)", "Petite family 2 (20 breakpoints/5 colonies)", "Petite family 3 (20 breakpoints/7 colonies)", "Petite family 4 (6 breakpoints/6 colonies)", "Petite family 5 (3 breakpoints/3 colonies)", "Petite family 6 (3 breakpoints/3 colonies)", "Petite family 7 (6 breakpoints/3 colonies)", "Petite family 8 (6 breakpoints/3 colonies)", "Petite family 9 (2 breakpoints/2 colonies)", "Grande colonies (3 breakpoints/10 colonies)"])

#
# #legend creation
# green_patch = mpatches.Patch(color='green', label='Negative strand origins')
# blue_patch = mpatches.Patch(color='blue', label='Positive strand origins')
# plt.legend(ncol=1, handles=[green_patch, blue_patch], loc='upper left', bbox_to_anchor=(1.0,0.5), fontsize=size)

figure = plt.gcf()
figure.set_size_inches(25,20)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir +os.sep

folder_name = os.sep + "junction_dot_plot" + os.sep
path = currDir + folder_name

print path
if not os.path.exists(path):
    os.makedirs(path)
    print path

#plt.savefig(path +  "normal_junctions", bbox_inches="tight", DPI=300)

#plt.savefig(path + "both_INREPEATS_junctions_with_families", bbox_inches="tight", dpi=300)
plt.savefig(path + "both_INREPEATS_junctions_with_families.svg", bbox_inches="tight")









