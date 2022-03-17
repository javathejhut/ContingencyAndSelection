import os
import sys
import glob
import pickle

import matplotlib.colors
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import random
import math
import scipy

# start by importing exp fit parameters
exp_fit_parameters_dict = {}

with open("readlength_exp_fit_parameters_petites.txt", 'r') as fparameters:
    for line in fparameters:
        name = line.strip().split()[0]
        location = float(line.strip().split()[1])
        scale = float(line.strip().split()[2])

        exp_fit_parameters_dict[name] = [location, scale]

#read in mt ref file
reference = ""
with open("chrM.fa", 'r') as ref_file:
    for line in ref_file:
        reference+= line.strip()

size = 50
# ori1, 2, 3, .., 8
origin_locations_forlines = [[4312, 4012], [32231, 32501], [54840, 54567], [56832, 56567], [82329, 82600],
                             [45227, 47927], [30220, 30594], [12780, 12510]]
origin_locations_active = [[4312, 4012], [32231, 32501], [54840, 54567],
                           [82329, 82600]]  # ori1, 2, 3, 5 "active" according to bernardi
mixed_samples = ["3I1", "6I1", "13I1", "17I1", "23I1", "8I1", "19I1"]
currDir = os.path.dirname(os.path.realpath(__file__))
folder_name = '/pikl_pipeline_output_petites/'
path = os.path.normpath(os.path.join(currDir, os.pardir)) + folder_name

sorted_primary_edges = {}
sorted_alternate_edges = {}

primary_repeat_lengths = {}
primary_repeat_lengths_double_mixed = {}

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


alternate_monomer_counts = {}
alternate_repeat_lengths = {}
primary_monomer_counts = {}

for file in os.listdir(path):
    with open(path + file, 'rb') as f:
        sample_name = file.split('_')[0]

        per_sample_pikl = pickle.load(f)

        primary_aln_key_count = per_sample_pikl[
            5]  # [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]
        alternate_aln_key_count = per_sample_pikl[
            6]  # [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]

        sorted_primary_edges[sample_name] = [sorted([ele[0], ele[1]]) for ele in primary_aln_key_count[1]]

        alternate_monomer_counts[sample_name] = []
        alternate_repeat_lengths[sample_name] = []
        primary_monomer_counts[sample_name] = []

        primary_monomer_counts[sample_name].append(primary_aln_key_count[-1])

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

        # compute primary repeat lengths from these
        repeat_length = 0
        for ele in sorted_primary_edges[sample_name]:
            repeat_length += abs(int(ele[1]) - int(ele[0]))

        primary_repeat_lengths[sample_name] = repeat_length

        #compute alternate edges
        sorted_alternate_edges[sample_name] = []
        for ele in alternate_aln_key_count:
            if ele[-1] != 0:
                structure_edges = []
                for aln in ele[1]:
                    to_append = sorted([aln[0], aln[1]])
                    structure_edges.append(to_append)
                sorted_alternate_edges[sample_name].append(structure_edges)
                # for alternate junctions (modified to be any monomer count)

        for ele in alternate_aln_key_count:
            # alternate_monomer_counts = []
            if ele[-1] != 0:
                alternate_monomer_counts[sample_name].append(ele[-1])

        for ele in alternate_aln_key_count:
            if ele[-1] != 0:
                repeat_lengths = 0
                for junc in ele[1]:
                    repeat_lengths += abs(junc[0] - junc[1])

                alternate_repeat_lengths[sample_name].append(repeat_lengths)

        # print sample_name, alternate_aln_key_count, sorted_alternate_edges[sample_name]
        # print sample_name, primary_aln_key_count, sorted_primary_edges[sample_name]

#computing alternate mt fractions (normalized)
def normalizer(some_mu, some_beta, some_L):
    return math.exp(-1 * (some_L - some_mu) / some_beta) + ((some_beta + some_mu) - (some_L + some_beta
                                                                                     ) * math.exp(
        -(some_L - some_mu) / some_beta)) / some_L

# now normalize counts for Petites (monomer counting instead)
petite_alternate_mt_fractions_mc = {}

for key in sorted_primary_edges.keys():
    # print key
    mu = float(exp_fit_parameters_dict[key][0])  # location of exponential
    beta = float(exp_fit_parameters_dict[key][1])  # scale of exponential

    L = primary_repeat_lengths[key]

    normalized_primary_content = primary_monomer_counts[key][0] * L / normalizer(mu, beta, L)

    alternate_normalized_content = []
    for i in range(len(alternate_monomer_counts[key])):
        L = alternate_repeat_lengths[key][i]

        normalized_alternate_content = alternate_monomer_counts[key][i] * L / normalizer(mu, beta, L)

        alternate_normalized_content.append(normalized_alternate_content)

    petite_alternate_mt_fractions_mc[key] = list(np.array(alternate_normalized_content)/(float(sum(alternate_normalized_content)) + normalized_primary_content))


# for alternates, return origin density
origin_fraction_primary = {}
origin_fraction_alternate = {}

edges_fraction_alternate = {}

origin_fraction_primary_active = {}
origin_fraction_alternate_active = {}

origin_fraction_wt_active = return_fractional_overlap([[0, 85779]], origin_locations_active)
origin_fraction_wt = return_fractional_overlap([[0, 85779]], origin_locations_forlines)

for key in sorted_alternate_edges:
    if len(sorted_alternate_edges[key]) > 0:

        origin_fraction_primary[key] = [return_fractional_overlap(sorted_primary_edges[key], origin_locations_forlines)]
        origin_fraction_alternate[key] = []
        edges_fraction_alternate[key] = []

        origin_fraction_primary_active[key] = [
            return_fractional_overlap(sorted_primary_edges[key], origin_locations_active)]
        origin_fraction_alternate_active[key] = []

        for structure in sorted_alternate_edges[key]:
            origin_fraction_alternate[key].append(return_fractional_overlap(structure, origin_locations_forlines))
            origin_fraction_alternate_active[key].append(return_fractional_overlap(structure, origin_locations_active))
            edges_fraction_alternate[key].append(structure)

        # print key, origin_fraction_primary[key], origin_fraction_alternate[key], origin_fraction_wt
        # print key, origin_fraction_primary_active[key], origin_fraction_alternate_active[key], origin_fraction_wt_active

#print origin_fraction_alternate
to_investigate = []
to_investigate_names = []
to_investigate_frequency = []
for key in origin_fraction_alternate.keys():
    if len(origin_fraction_alternate[key])>0:
        for i, fraction in enumerate(origin_fraction_alternate[key]):
            if fraction ==0:
                to_investigate.append(edges_fraction_alternate[key][i])
                to_investigate_names.append(key)
                to_investigate_frequency.append(petite_alternate_mt_fractions_mc[key][i])
                #print key, fraction, edges_fraction_alternate[key][i]

mean_w = 10
fraction_range = np.linspace(0, 1.0, mean_w +1)

#generate reference GC window content
ref_GC_content = []
ref_GC_fraction = []
for i in range(0, len(reference) - mean_w + 1):
    in_origin=False
    for loc in origin_locations_forlines:
        if i in range(loc[0], loc[1]) or i+mean_w in range(loc[0], loc[1]):
            in_origin=True

    #if not in_origin:
    window = reference[i: i + mean_w]

    avg_GC_ref = float(window.count('C') + window.count('G')) / len(window)
    ref_GC_content.append(avg_GC_ref)

    # print sim_GC_content
for frac in fraction_range:
    # print frac, sim_GC_content,
    ref_GC_fraction.append(len([ele for ele in ref_GC_content if ele >= frac]) / float(len(ref_GC_content)))

#now do for sample :)
labeled=False
index = 0

cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(min(to_investigate_frequency), max(to_investigate_frequency))
num_lines_above_red = [0]*len(fraction_range)

for sample in to_investigate:

    name = to_investigate_names[index]
    #print name, sample, petite_alternate_mt_fractions_mc[name], to_investigate_frequency[index]


    test_sample = reference[int(sample[0][0])-1:int(sample[0][1])+1]

    GC_content = []
    frac_gc_content_data = []
    for i in range(0, len(test_sample) - mean_w + 1):
        window = test_sample[i: i + mean_w]

        avg_GC = float(window.count('C') + window.count('G')) / len(window)
        GC_content.append(avg_GC)

    for frac in fraction_range:
        # print frac, sim_GC_content,
        frac_gc_content_data.append(len([ele for ele in GC_content if ele >= frac]) / float(len(GC_content)))
        #print frac
        if frac ==0.8 and frac_gc_content_data[-1] < ref_GC_fraction[-3]:
            print name, frac_gc_content_data[-1], "special"

        elif frac==0.8:
            print name, frac_gc_content_data[-1]

    for i, frac in enumerate(fraction_range):
        if frac_gc_content_data[i] > ref_GC_fraction[i]:
            num_lines_above_red[i]+=1

plt.plot(fraction_range, np.array(num_lines_above_red)/9.0, 'g-', linewidth=4)

    # if labeled==False:
    #     #plt.plot(fraction_range, frac_gc_content_data, color=cmap(norm(to_investigate_frequency[index])), linewidth=3,  label="GC content in alternate structures w/o origins")
    #     plt.plot(fraction_range, frac_gc_content_data, 'k', linewidth=4,
    #              label="GC content in alternate structures w/o origins")
    # else:
    #     #plt.plot(fraction_range, frac_gc_content_data,  color=cmap(norm(to_investigate_frequency[index])), linewidth=3 )
    #     plt.plot(fraction_range, frac_gc_content_data, 'k', linewidth=4)
    #
    # labeled = True
    # index += 1
    #plt.xlim(0.6, 1)
    #plt.plot(range(len(GC_content)), GC_content)
    #plt.show()
fill_negative = plt.fill_between(np.arange(0.6, 1.1, 0.01), -0.1, 8. , color='gray', alpha= 0.3)
#plt.plot([0.6]*81, np.arange(0,8.1,0.1 ), 'b--', linewidth=4, alpha= 0.5, label='Window GC fraction that defines GC clusters')
#plt.plot(fraction_range, ref_GC_fraction, 'r--', linewidth=4, label="Average GC content across mt reference")
#plt.xlabel("GC content fraction in %i bp window"%mean_w, fontsize=size)
plt.xticks([0.0, 0.2,  0.4, 0.6,  0.8,  1.0], [0.0, 0.2,  0.4, 0.6,  0.8,  1.0],fontsize=size, rotation=0)


plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=size, rotation=0)
plt.xlim(-0.05,1)
plt.ylim(0, 1)

ax = plt.gca()
ax.tick_params(axis='x', pad=20)

#plt.legend(fontsize=size-8, loc="upper right")
plt.ylabel("Fraction of structures\n enriched in GC", fontsize=size)

figure = plt.gcf()
figure.set_size_inches(14, 8)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep
#
# folder_name = os.sep + "ori_distance_dists" + os.sep
# path = currDir + folder_name
#
# if not os.path.exists(path):
#     os.makedirs(path)
#
# #save histogram figure
# #plt.show()
plt.savefig(path + "zero_ori_alternate_GC_coloured_by_freq_inset.png", bbox_inches="tight", dpi=300, transparent=False, linewidth=2)
#
# figure.clf()







        #print  petite_alternate_mt_fractions_mc[key], origin_fraction_alternate[key]

# for key in origin_fraction_alternate.keys():
#     if len(origin_fraction_alternate[key])>1:
#         print  petite_alternate_mt_fractions_mc[key], origin_fraction_alternate[key]
#
# fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
# #fig.suptitle('Sharing x per column, y per row')
#
# axes_to_plot = [ax1, ax2, ax3, ax4, ax5, ax6]
# axes_index = 0
# for key in origin_fraction_alternate.keys():
#     if len(origin_fraction_alternate[key])>1:
#
#         axes_to_plot[axes_index].scatter(origin_fraction_alternate[key], petite_alternate_mt_fractions_mc[key])
#         #z = np.polyfit(x=origin_fraction_alternate[key], y= petite_alternate_mt_fractions_mc[key], deg=1)
#
#         #xx = np.linspace(0, 1, 40)
#         pearson_r = scipy.stats.pearsonr(origin_fraction_alternate[key], petite_alternate_mt_fractions_mc[key])
#         axes_to_plot[axes_index].text(0.5, 0.5, "$R^2$ = %2.2f"%pearson_r[0], size=size)
#         #axes_to_plot[axes_index].plot(xx, z[0]*xx + z[1]*xx, 'k')
#         axes_index += 1
#
# for ax in fig.get_axes():
#     ax.label_outer()
#     ax.set_xlim(0, 1)
#     ax.set_ylim(0, 1)
#     ax.tick_params(labelsize=size)
#
# plt.xticks(fontsize=size, rotation=0)
# plt.yticks(fontsize=size, rotation=0)
#
# fig.add_subplot(111, frameon=False)
# # hide tick and tick label of the big axes
# plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
# plt.xlabel("Origin content fraction \nof alternate structures", fontsize=size, labelpad=20)
# plt.ylabel("Relative frequency of \nalternate structures", fontsize=size, labelpad=20)
# #plt.show()
#
# figure = plt.gcf()
# figure.set_size_inches(10, 12)
#
# currDir = os.path.dirname(os.path.realpath(__file__))
# path = currDir + os.sep
#
# plt.savefig(path + "origin_fraction_freq_correlation.png", bbox_inches="tight", dpi=300, format='PNG')
#
#
# plt.clf()
#
# plot_alternate = origin_fraction_alternate
# plot_primary = origin_fraction_primary
# plot_wt_line = origin_fraction_wt
#
# # do some sorting please, to make this look reasonable
#
# # converting to list
#
# # plot_primary, plot_alternate = (list(t) for t in zip(*sorted(zip(plot_primary, plot_alternate), reverse=False)))
#
# # print plot_primary
# # now do some fancy plotting, also add line for WT version
# labels_primary = []
# labels_alternate = []
#
# hues_primary = []
# hues_alternate = []
#
# for i, key in enumerate(plot_primary):
#     labels_primary.extend([key] * len(plot_primary[key]))
#     labels_alternate.extend([key] * len(plot_alternate[key]))
#
#     hues_primary.extend(["primary"] * len(plot_primary[key]))
#     hues_alternate.extend(["alternate"] * len(plot_alternate[key]))
#
# # print labels
#
# size = 25
# toplot_primary = []
# toplot_alternate = []
#
# for key in plot_primary:
#     toplot_primary.extend(plot_primary[key])
#     toplot_alternate.extend(plot_alternate[key])
#
# # plotting line
# plt.plot(range(0, 16), [plot_wt_line] * 16, 'g--', linewidth=2)
#
# for i in range(0, 16):
#     plt.plot([i] * len(np.arange(-0.1, 1.1, 0.05)), np.arange(-0.1, 1.1, 0.05), '-', color='gray', linewidth=2,
#              alpha=0.5)
# plt.plot(range(-1, 17), [plot_wt_line] * 18, 'g--', linewidth=2)
#
# d_primary = {'Sample': labels_primary, 'Fractional origin content per structure': toplot_primary,
#              "Structure type": hues_primary}
# d_alternate = {'Sample': labels_alternate, 'Fractional origin content per structure': toplot_alternate,
#                "Structure type": hues_alternate}
#
# df_primary = pd.DataFrame(d_primary).sort_values('Fractional origin content per structure')
# df_alternate = pd.DataFrame(d_alternate)
#
# df = pd.concat([df_primary, df_alternate])
#
# # sns.violinplot(data = df, x='Junction', y='Distance to closest origin (bp)', color="0.8", cut=0, inner=None)
# ax = sns.stripplot(data=df, x='Sample', y='Fractional origin content per structure', jitter=True, size=10,
#                    hue="Structure type")
# # ax = sns.boxplot(data=df, x='Excision rate (per generation)', y='Alternate mt content fraction',palette="Set2",fliersize=0, hue="Petite lineages")
# # plt.scatter(x=range(len(means)), y=means, c='k')
# # plt.ylim(0.0, 0.025)
# handles, legend_labels = ax.get_legend_handles_labels()
# # print legend_labels
# # plt.xticks([0, 1,2,3,4,5,6,7,8,9],[0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05], fontsize=size)
# legend = plt.legend(handles[-2:], legend_labels[-2:], fontsize=size, loc="upper left")
# plt.xticks(fontsize=size, rotation=270)
# plt.yticks(fontsize=size, rotation=0)
# plt.ylim(-0.05, 1.0)
# plt.xlabel("Sample", fontsize=size)
# plt.ylabel("Origin content fraction", fontsize=size)
# plt.setp(plt.gca().get_legend().get_texts(), fontsize=size)
# legend.set_title("Structure Type", prop={'size': size})
# figure = plt.gcf()
# figure.set_size_inches(12, 10)
#
# currDir = os.path.dirname(os.path.realpath(__file__))
# path = currDir + os.sep
#
# plt.savefig(path + "origin_fraction_primary_alternate.png", bbox_inches="tight", dpi=300, format='PNG')

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
