import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import random
import math
import scipy

'''Gather all junctions/types in pickled output across all strains'''

# start by importing exp fit parameters
exp_fit_parameters_dict = {}

with open("readlength_exp_fit_parameters_petites.txt", 'r') as fparameters:
    for line in fparameters:
        name = line.strip().split()[0]
        location = float(line.strip().split()[1])
        scale = float(line.strip().split()[2])

        exp_fit_parameters_dict[name] = [location, scale]

# construct all dicts... for all types
alternate_repeat_lengths = {}
alternate_monomer_counts = {}

primary_repeat_lengths = {}
primary_monomer_counts = {}

print(os.path.normpath(os.path.join(os.getcwd(), os.pardir)))
# for petites
for sample_pikl_file in glob.glob(os.path.normpath(os.path.join(os.getcwd(), os.pardir))+'/pikl_pipeline_output_petites/*.pikl'):
    pikl_name = sample_pikl_file.strip().split('/')[-1].split('_')[0]

    with open(sample_pikl_file, 'rb') as fopen:
        per_sample_pikl = pickle.load(fopen)

        mixed_alignment_count_dict = per_sample_pikl[
            0]  # key: LH format, value: [[aln start, aln end, aln start std, aln end std], count]
        example_mixed_read_dict = per_sample_pikl[
            1]  # key: LH format, value: [[[LH representation]], [[alphabet representation]], [[aln1, aln end 1, aln1 std, aln end 1 std], [aln2, aln end2, aln2 std, aln 2 end std]]]
        primary_trans_dict = per_sample_pikl[
            2]  # key: (transition tuple), value: [count no sign transition, count sign transition]

        inferred_primary = per_sample_pikl[
            3]  # [[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]
        primary_lcs_suggestion_dict = per_sample_pikl[
            4]  # key: (lH repeat format), value : [[orientation rep], [alphabet rep], [[aln1, aln1 end, aln1 std, aln1 end std], [aln2 , aln2 end, aln2 std, aln2 end std]]]
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

        '''do list initialization'''

        alternate_monomer_counts[pikl_name] = []
        alternate_repeat_lengths[pikl_name] = []

        primary_monomer_counts[pikl_name] = []
        primary_repeat_lengths[pikl_name] = []

        # for primary junctions
        primary_junction_counts = []
        primary_junctions = []
        for junc in primary_aln_key_count[0]:
            primary_junction_counts.append(junction_counts_dict[junc])
            primary_junctions.append(junc)

        primary_monomer_counts[pikl_name].append(primary_aln_key_count[-1])

        repeat_lengths = 0

        if not pikl_name == '10I2' and not pikl_name == "11I2" and not pikl_name == "20I2" and not pikl_name == "21I2":
            for ele in primary_aln_key_count[1]:
                repeat_lengths += abs(float(ele[0]) - float(ele[1]))
        else:
            for ele in primary_aln_key_count[1]:
                repeat_lengths += 85779 - (abs(float(ele[0]) - float(ele[1])))

        primary_repeat_lengths[pikl_name].append(repeat_lengths)

        # for alternate junctions (modified to be any monomer count)
        for ele in alternate_aln_key_count:
            # alternate_monomer_counts = []
            if ele[-1] != 0:
                alternate_monomer_counts[pikl_name].append(ele[-1])

        for ele in alternate_aln_key_count:
            if ele[-1] != 0:
                repeat_lengths = 0
                for junc in ele[1]:
                    repeat_lengths += abs(junc[0] - junc[1])

                alternate_repeat_lengths[pikl_name].append(repeat_lengths)

        # adjustments for 12I1 outlier
        if pikl_name == "12I1":
            alternate_monomer_counts[pikl_name].append(primary_aln_key_count[-1])
            alternate_repeat_lengths[pikl_name].append(primary_repeat_lengths[pikl_name][0])

# print alternate_monomer_counts
# print alternate_repeat_lengths

def normalizer(some_mu, some_beta, some_L):
    return math.exp(-1 * (some_L - some_mu) / some_beta) + ((some_beta + some_mu) - (some_L + some_beta
                                                                                     ) * math.exp(
        -(some_L - some_mu) / some_beta)) / some_L

# now normalize counts for Petites (monomer counting instead)
petite_alternate_mt_fractions_mc = {}
petite_mt_fractions_mc = {}

for key in primary_monomer_counts.keys():
    # print key
    mu = float(exp_fit_parameters_dict[key][0])  # location of exponential
    beta = float(exp_fit_parameters_dict[key][1])  # scale of exponential

    alternate_normalized_content = []
    normalized_primary_content = 4*50000/normalizer(mu, beta, 50000)
    print primary_monomer_counts[key]
    for i in range(len(alternate_monomer_counts[key])):
        L = alternate_repeat_lengths[key][i]

        normalized_alternate_content = alternate_monomer_counts[key][i] * L / normalizer(mu, beta, L)

        alternate_normalized_content.append(normalized_alternate_content)

    petite_alternate_mt_fractions_mc[key] = list(np.array(alternate_normalized_content)/float(sum(alternate_normalized_content)))

    print key, list(np.array(alternate_normalized_content)/float(sum(alternate_normalized_content) + normalized_primary_content)), normalized_primary_content/float(sum(alternate_normalized_content) + normalized_primary_content)

    #petite_mt_fractions_mc[key] =
    # petite_alternate_mt_fractions_mc[key] = sum(alternate_normalized_content) / float(
    #     sum(alternate_normalized_content) + normalized_primary_content)

#print petite_alternate_mt_fractions_mc

size = 25
# ori1, 2, 3, .., 8
origin_locations_forlines = [[4312, 4012], [32231, 32501], [54840, 54567], [56832, 56567], [82329, 82600],
                             [45227, 47927], [30220, 30594], [12780, 12510]]
origin_locations_active = [[4312, 4012], [32231, 32501], [54840, 54567],
                           [82329, 82600]]  # ori1, 2, 3, 5 "active" according to bernardi
mixed_samples = ["3I1", "4I1", "5I1","6I1", "10I1", "13I1", "17I1", "18I1", "20I1" , "23I1", "5I2",
  "8I1", "19I1"]
currDir = os.path.dirname(os.path.realpath(__file__))
folder_name = '/pikl_pipeline_output_petites/'
path = os.path.normpath(os.path.join(currDir, os.pardir)) + folder_name
sorted_primary_edges = {}
sorted_alternate_edges = {}

def return_fractional_overlap(query_list, reference_list):
    query_length = 0
    # perform pairwise query to reference alignment test
    total_query_overlap_bp = 0
    # print query_list
    for q_list in query_list:
        per_q_overlap = []
        query_length += abs(min(q_list[0], q_list[1]) - max(q_list[0], q_list[1]))
        for r_list in reference_list:
            q_range = set(range(min(q_list[0], q_list[1]), max(q_list[0], q_list[1]) + 1))
            r_range = set(range(min(r_list[0], r_list[1]), max(r_list[0], r_list[1]) + 1))

            # symm_diff = (q_range ^ r_range) & q_range
            overlap = (q_range & r_range)
            per_q_overlap.append(len(overlap))

        total_query_overlap_bp += sum(per_q_overlap)
    # print total_query_overlap_bp, query_length
    return float(total_query_overlap_bp) / query_length

for file in os.listdir(path):
    with open(path+file, 'rb') as f:
        sample_name = file.split('_')[0]

        per_sample_pikl = pickle.load(f)

        primary_aln_key_count = per_sample_pikl[
            5]  # [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]
        alternate_aln_key_count = per_sample_pikl[
            6]  # [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]

        sorted_primary_edges[sample_name] = [sorted([ele[0], ele[1]]) for ele in primary_aln_key_count[1]]

        '''Here we want edges that define alignments in ONE repeat unit'''
        if sample_name == '12I1':
            sorted_primary_edges[sample_name] = [[sorted_primary_edges[sample_name][1][0],
                                                  sorted_primary_edges[sample_name][0][0]]]

        if sample_name == "10I2" or sample_name == "11I2" or sample_name == "20I2" or sample_name == "21I2":
            # print sorted_primary_edges[sample_name]
            sorted_primary_edges[sample_name] = [[0, sorted_primary_edges[sample_name][0][0]],
                                                 [sorted_primary_edges[sample_name][0][1], 85779]]

        # handle mixed structure cases, we want alignments for just longest and smallest alignments
        if sample_name in mixed_samples:
            sizes = [abs(ele[1] - ele[0]) for ele in sorted_primary_edges[sample_name]]
            max_index = sizes.index(max(sizes))
            min_index = sizes.index(min(sizes))

            sorted_primary_edges[sample_name] = [sorted_primary_edges[sample_name][min_index],
                                                 sorted_primary_edges[sample_name][max_index]]

        sorted_alternate_edges[sample_name] = []
        for ele in alternate_aln_key_count:
            if ele[-1] != 0:
                structure_edges = []
                for aln in ele[1]:
                    to_append = sorted([aln[0], aln[1]])
                    structure_edges.append(to_append)
                sorted_alternate_edges[sample_name].append(structure_edges)

        # print sample_name, alternate_aln_key_count, sorted_alternate_edges[sample_name]
        # print sample_name, primary_aln_key_count, sorted_primary_edges[sample_name]


# for alternates, return origin density
origin_fraction_primary = {}
origin_fraction_alternate = {}

origin_fraction_primary_active = {}
origin_fraction_alternate_active = {}

origin_fraction_wt_active = return_fractional_overlap([[0, 85779]], origin_locations_active)
origin_fraction_wt = return_fractional_overlap([[0, 85779]], origin_locations_forlines)

for key in sorted_alternate_edges:
    if len(sorted_alternate_edges[key]) > 0:

        origin_fraction_primary[key] = [return_fractional_overlap(sorted_primary_edges[key], origin_locations_forlines)]
        origin_fraction_alternate[key] = []

        origin_fraction_primary_active[key] = [
            return_fractional_overlap(sorted_primary_edges[key], origin_locations_active)]
        origin_fraction_alternate_active[key] = []

        for structure in sorted_alternate_edges[key]:
            origin_fraction_alternate[key].append(return_fractional_overlap(structure, origin_locations_forlines))
            origin_fraction_alternate_active[key].append(return_fractional_overlap(structure, origin_locations_active))

        # print key, origin_fraction_primary[key], origin_fraction_alternate[key], origin_fraction_wt
        # print key, origin_fraction_primary_active[key], origin_fraction_alternate_active[key], origin_fraction_wt_active

#print origin_fraction_alternate

for key in petite_alternate_mt_fractions_mc.keys():
    if len(petite_alternate_mt_fractions_mc[key])>1:
        print  petite_alternate_mt_fractions_mc[key], origin_fraction_alternate[key]

#fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
#fig.suptitle('Sharing x per column, y per row')

axes_to_plot = [ax1, ax2, ax3, ax4, ax5, ax6]
text_loc = [[5000, 0.4], [5000, 0.15], [6000, 0.3], [700, 0.51], [1500, 0.6], [6000, 0.5]]
axes_index = 0
for key in origin_fraction_alternate.keys():
    if len(origin_fraction_alternate[key])>1:

        #axes_to_plot[axes_index].scatter(origin_fraction_alternate[key], petite_alternate_mt_fractions_mc[key])
        axes_to_plot[axes_index].scatter(alternate_repeat_lengths[key], petite_alternate_mt_fractions_mc[key], s=50, color='k')
        #z = np.polyfit(x=origin_fraction_alternate[key], y= petite_alternate_mt_fractions_mc[key], deg=1)

        #xx = np.linspace(0, 1, 40)
        #pearson_r = scipy.stats.pearsonr(origin_fraction_alternate[key], petite_alternate_mt_fractions_mc[key])
        pearson_r = scipy.stats.pearsonr(alternate_repeat_lengths[key], petite_alternate_mt_fractions_mc[key])
        #axes_to_plot[axes_index].text(0.5, 0.5, "$R^2$ = %2.2f"%pearson_r[0], size=size)
        axes_to_plot[axes_index].text(text_loc[axes_index][0], text_loc[axes_index][1], "$R^2$ = %2.2f" % pearson_r[0], size=size)
        #axes_to_plot[axes_index].plot(xx, z[0]*xx + z[1]*xx, 'k')
        axes_index += 1

for ax in fig.get_axes():
    #ax.label_outer()
    #ax.set_xlim(0, 1)
    #ax.set_ylim(0, 1)
    ax.tick_params(labelsize=size)

# for i, ax in enumerate(axes_to_plot):
#     if i!=0 and i!=4:
#         ax.tick_params(labelleft=False)

plt.xticks(fontsize=size, rotation=0)
plt.yticks(fontsize=size, rotation=0)

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel("Repeat unit length (bp)", fontsize=size, labelpad=20)
plt.ylabel("Relative frequency of alternate structures\nwithin petite strains", fontsize=size, labelpad=38)
#plt.show()

figure = plt.gcf()
figure.set_size_inches(21, 10)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

plt.savefig(path + "origin_fraction_freq_correlation.png", bbox_inches="tight", dpi=300, format='PNG')


plt.clf()

plot_alternate = origin_fraction_alternate
plot_primary = origin_fraction_primary
plot_wt_line = origin_fraction_wt

# do some sorting please, to make this look reasonable

# converting to list

# plot_primary, plot_alternate = (list(t) for t in zip(*sorted(zip(plot_primary, plot_alternate), reverse=False)))

# print plot_primary
# now do some fancy plotting, also add line for WT version
labels_primary = []
labels_alternate = []

hues_primary = []
hues_alternate = []

for i, key in enumerate(plot_primary):
    labels_primary.extend([key] * len(plot_primary[key]))
    labels_alternate.extend([key] * len(plot_alternate[key]))

    hues_primary.extend(["primary"] * len(plot_primary[key]))
    hues_alternate.extend(["alternate"] * len(plot_alternate[key]))

# print labels

size = 32
toplot_primary = []
toplot_alternate = []

for key in plot_primary:
    toplot_primary.extend(plot_primary[key])
    toplot_alternate.extend(plot_alternate[key])

# plotting line
plt.plot(range(0, 16), [plot_wt_line] * 16, 'g--', linewidth=2)

for i in range(0, 16):
    plt.plot([i] * len(np.arange(-0.1, 1.1, 0.05)), np.arange(-0.1, 1.1, 0.05), '-', color='gray', linewidth=2,
             alpha=0.5)
plt.plot(range(-1, 17), [plot_wt_line] * 18, 'g--', linewidth=2, label='Wild-type mtDNA')

d_primary = {'Sample': labels_primary, 'Fractional origin content per structure': toplot_primary,
             "Structure type": hues_primary}
d_alternate = {'Sample': labels_alternate, 'Fractional origin content per structure': toplot_alternate,
               "Structure type": hues_alternate}

df_primary = pd.DataFrame(d_primary).sort_values('Fractional origin content per structure')
df_alternate = pd.DataFrame(d_alternate)

df = pd.concat([df_primary, df_alternate])

# sns.violinplot(data = df, x='Junction', y='Distance to closest origin (bp)', color="0.8", cut=0, inner=None)
ax = sns.stripplot(data=df, x='Sample', y='Fractional origin content per structure', jitter=True, size=10,
                   hue="Structure type")
# ax = sns.boxplot(data=df, x='Excision rate (per generation)', y='Alternate mt content fraction',palette="Set2",fliersize=0, hue="Petite lineages")
# plt.scatter(x=range(len(means)), y=means, c='k')
# plt.ylim(0.0, 0.025)
handles, legend_labels = ax.get_legend_handles_labels()
print legend_labels
# print legend_labels
# plt.xticks([0, 1,2,3,4,5,6,7,8,9],[0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05], fontsize=size)
#legend = plt.legend([handles[0]] + handles[-2:], [legend_labels[0]] + legend_labels[-2:], fontsize=size, loc="upper left")
legend = plt.legend([handles[0]] + handles[-2:], ['Wild-type mtDNA', 'Primary Petite structure', 'Alternate Petite structure'], fontsize=size, loc="upper left")
plt.xticks([0,1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],fontsize=size-5)
plt.yticks(fontsize=size, rotation=0)
plt.ylim(-0.05, 1.05)
plt.xlabel("Petite colony", fontsize=size)
plt.ylabel("Origin content fraction", fontsize=size)
plt.setp(plt.gca().get_legend().get_texts(), fontsize=size)
#legend.set_title("Structure Type", prop={'size': size})
figure = plt.gcf()
figure.set_size_inches(12, 10)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

plt.savefig(path + "origin_fraction_primary_alternate.png", bbox_inches="tight", dpi=300, format='PNG')

# iterate through alternate junctions, compare all to all primary of same key
# fractional_distances_to_plot = []
# for key_alternate in sorted_alternate_edges.keys():
#
#     for key_primary in sorted_primary_edges.keys():
#         if key_alternate == key_primary:
#
#             alternate_edges = sorted_alternate_edges[key_alternate]
#             primary_edges = sorted_primary_edges[key_alternate]
#
#             #for every alternate, compute smallest fractional distance to edge
#             for edge_a in alternate_edges:
#                 frac_distances = []
#                 for edge_p in primary_edges:
#
#                     #two samples where complement is actually primary structure
#                     if key_alternate == "10I2" or key_alternate == "11I2":
#                         frac_dist = min(abs(edge_a[0]-edge_p[1])/float(85779- abs(edge_p[1]-edge_p[0])),
#                                         abs(edge_a[1]-85779)/float(85779-abs(edge_p[1]-edge_p[0])))
#                     else:
#                         frac_dist = min(abs(edge_a[0] - edge_p[0]) / float(abs(edge_p[1] - edge_p[0])),
#                                         abs(edge_a[1] - edge_p[1]) / float(abs(edge_p[1] - edge_p[0])))
#                     frac_distances.append(frac_dist)
#
#                 if key_alternate=="12I1":
#                     print min(frac_distances)
#                 fractional_distances_to_plot.append(min(frac_distances))
#
# #print fractional_distances_to_plot
# #sns.stripplot(y = distances_to_plot, jitter=True, orient='h')
# sns.swarmplot(y = fractional_distances_to_plot,  orient='h', size=8, color='b')
#
# #n, bins, edges = plt.hist(distances_to_plot, bins= 100, density=True, alpha=0.5, label="Observed, mean=%2.2f" %(np.median(distances_to_plot)))
#
# #n_alternate, bins_alternate, edges_alternate = plt.hist(distances_to_plot_alternate, bins= 20, density=True, alpha=0.5, label="Observed, mean=%2.2f" %(np.mean(distances_to_plot_alternate)))
#
# plt.title("Regions near excision sites involved in primary structures appear to be \n preferred in forming subsequent petite structures", fontsize=size)
# plt.xlabel("Fractional distance between primary alignment edge and closest alternate alignment edge", fontsize=size)
# #plt.text(1000, 0.0015, 'Observed: 66% of edges are within 580bp of an origin', fontsize=size)
# plt.text(0.05, 0.3, 'mean= %2.2f \nmedian= %2.2f '%(np.mean(fractional_distances_to_plot), np.median(fractional_distances_to_plot)), fontsize=size)
# #plt.ylabel("Frequency", fontsize=size)
# #plt.xlabel("Source of junctions")
# plt.xticks(fontsize=size, rotation=0)
# #plt.yticks(fontsize=size, rotation=0)
# #plt.xlim(0, 5000)
# #plt.legend(loc='upper right', fontsize=size)
# figure = plt.gcf()
# figure.set_size_inches(12, 10)
#
# currDir = os.path.dirname(os.path.realpath(__file__))
# path = currDir + os.sep
#
# folder_name = os.sep + "ori_distance_dists" + os.sep
# path = currDir + folder_name
#
# if not os.path.exists(path):
#     os.makedirs(path)
#
# #save histogram figure
# #plt.show()
# plt.savefig(path + "primary_alternate_edge_distances_relative.png", bbox_inches="tight", DPI=700)
#
# figure.clf()
