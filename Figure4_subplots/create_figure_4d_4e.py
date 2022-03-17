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
import matplotlib as mpl

def is_overlapped(a_start, a_end, b_start, b_end):
    query1_range = [int(a_start), int(a_end)]
    query2_range = [int(b_start), int(b_end)]
    #if partially overlapping
    if (int(b_start) in range(min(query1_range), max(query1_range))) \
            or (int(b_end) in range(min(query1_range), max(query1_range))):
        return True
    #if b is encompassing a
    elif (int(a_start) in range(min(query2_range), max(query2_range))) \
            or (int(a_end) in range(min(query2_range), max(query2_range))):
        return True
    else:
        return False

'''Gather all junctions/types in pickled output across all strains'''
size=32

currDir = os.path.dirname(os.path.realpath(__file__))
folder_name = '/pikl_pipeline_output_petites/'
path = os.path.normpath(os.path.join(currDir, os.pardir)) + folder_name

sorted_primary_edges = {}
sorted_primary_edges_std = {} #standard deviation in edge locations

sorted_alternate_edges = {}
sorted_alternate_edges_std = {} #standard deviation in alternage edge locations

for file in os.listdir(path):
    with open(path + file, 'rb') as f:
        sample_name = file.split('_')[0]

        per_sample_pikl = pickle.load(f)

        primary_aln_key_count = per_sample_pikl[5]  # [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]
        alternate_aln_key_count = per_sample_pikl[6]  # [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]

        sorted_primary_edges[sample_name] = [sorted([ele[0],ele[1]]) for ele in primary_aln_key_count[1]]
        sorted_primary_edges_std[sample_name] = [sorted([ele[2], ele[3]]) for ele in primary_aln_key_count[1]]

        if sample_name == '12I1':

            sorted_primary_edges[sample_name] = [[sorted_primary_edges[sample_name][1][0],
                                                       sorted_primary_edges[sample_name][0][0]]]
            sorted_primary_edges_std[sample_name] = [[1000,1000]]

            print sorted_primary_edges[sample_name]
            alternate_aln_key_count.append(primary_aln_key_count)
            #print alternate_aln_key_count
            
        sorted_alternate_edges[sample_name] = []
        sorted_alternate_edges_std[sample_name] = []

        for ele in alternate_aln_key_count:

            #if we have counts for this structure
            if ele[-1]>0:
                edges_for_one_structure = []
                edges_for_one_structure_std = []
                for aln in ele[1]:
                    to_append = [aln[0], aln[1]]
                    to_append_std = [aln[2], aln[3]]
                    edges_for_one_structure.extend(to_append)
                    edges_for_one_structure_std.extend(to_append_std)

                sorted_alternate_edges[sample_name].append(edges_for_one_structure)
                sorted_alternate_edges_std[sample_name].append(edges_for_one_structure_std)

#iterate through alternate junctions, compare all to all primary of same key, and classify

fractional_distances_across_primary = []
fractional_distances_shared_primary = []
fractional_distances_inside_primary = []

distances_across_primary = []
distances_shared_primary = []
distances_inside_primary = []

for key_alternate in sorted_alternate_edges.keys():

    for key_primary in sorted_primary_edges.keys():
        if key_alternate == key_primary:

            alternate_edges = sorted_alternate_edges[key_alternate]
            alternate_edge_std = sorted_alternate_edges_std[key_alternate]

            primary_edges = sorted_primary_edges[key_alternate]
            primary_edges_std = sorted_primary_edges_std[key_alternate]

            #for every alternate, compute smallest fractional distance to edge
            for repeat_unit, repeat_unit_std in zip(alternate_edges, alternate_edge_std):
                frac_distances = []
                distances = []
                includes_primary=False
                for edge_a, edge_a_std in zip(repeat_unit, repeat_unit_std):
                    for edge_p, edge_p_std in zip(primary_edges, primary_edges_std):
                        #print edge_a, edge_a_std
                        if is_overlapped(edge_a - edge_a_std, edge_a + edge_a_std, edge_p[0] - edge_p_std[0], edge_p[0] + edge_p_std[0])\
                                or is_overlapped(edge_a - edge_a_std, edge_a + edge_a_std, edge_p[1] - edge_p_std[1], edge_p[1] + edge_p_std[1]):
                            includes_primary=True

                        #two samples where complement is actually primary structure
                        if key_alternate == "10I2" or key_alternate == "11I2" or key_alternate=="20I2" or key_alternate=="21I2":
                            frac_dist = min(abs(edge_a-edge_p[1])/float(85779- abs(edge_p[1]-edge_p[0])),
                                            abs(edge_a-85779)/float(85779-abs(edge_p[1]-edge_p[0])))
                            dist = min(abs(edge_a-edge_p[1]), abs(edge_a-85779))
                        else:
                            frac_dist = min(abs(edge_a - edge_p[0]) / float(abs(edge_p[1] - edge_p[0])),
                                            abs(edge_a - edge_p[1]) / float(abs(edge_p[1] - edge_p[0])))

                            dist = min(abs(edge_a - edge_p[0]), abs(edge_a - edge_p[1]))

                        frac_distances.append(frac_dist)
                        distances.append(dist)
                    if key_alternate=="3I1":
                        print "**", min(frac_distances), includes_primary
                if len(repeat_unit)>2 and includes_primary:
                    fractional_distances_across_primary.append(min(frac_distances))
                    distances_across_primary.append(min(distances))
                elif (len(repeat_unit)==2 and includes_primary) or key_primary=="4I1" or key_primary=="10I2" or key_primary=="11I2": #handle cases where alternates outside primary alignment, which necessarily includes primary edges
                    print key_primary
                    fractional_distances_shared_primary.append(min(frac_distances))
                else:
                    fractional_distances_inside_primary.append(min(frac_distances))

print np.mean(distances_across_primary)

print len(fractional_distances_across_primary)
print len(fractional_distances_shared_primary)
print len(fractional_distances_inside_primary)

print fractional_distances_inside_primary

#first create a bar plot showing the percentage of alternate structures in each class
total_count = float(len(fractional_distances_across_primary) + len(fractional_distances_shared_primary) + len(fractional_distances_inside_primary))
for_y_bar_fraction = [ len(fractional_distances_inside_primary), len(fractional_distances_across_primary) , len(fractional_distances_shared_primary)]
patches, texts, autotexts = plt.pie(for_y_bar_fraction, labels=["Type I", "Type II", "Type III"],  autopct=lambda(p): '{:.0f}%'.format(round(p)), colors=["#d62728","#9467bd","#17becf"], textprops={'fontsize': size}, startangle=90)

print total_count
print for_y_bar_fraction
for t in texts:
    t.set_fontsize(size)
    
#plt.ylabel("% of alternate structures\nwith this excision type", fontsize=size)
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.axis('equal')
figure = plt.gcf()
figure.set_size_inches(14, 9)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

folder_name = os.sep + "ori_distance_dists" + os.sep
path = currDir + folder_name

if not os.path.exists(path):
    os.makedirs(path)

#save histogram figure
#plt.show()
plt.savefig(path + "structure_classes_barchart_std_piechart.png", bbox_inches="tight", dpi=300)

figure.clf()

#pandas formatting for multi stripplot

shared_label = ["Type III" for i in range(len(fractional_distances_shared_primary))]
across_label = ["Type II" for i in range(len(fractional_distances_across_primary))]
within_label = ["Type I" for i in range(len(fractional_distances_inside_primary))]

labels = []
labels.extend(within_label)
labels.extend(across_label)
labels.extend(shared_label)

toplot = []
toplot.extend(fractional_distances_inside_primary)
toplot.extend(fractional_distances_across_primary)
toplot.extend(fractional_distances_shared_primary)

d = {'Class': labels, '% of alternate structures\nwith this excision type' : toplot}
df = pd.DataFrame(d)

palette = sns.color_palette(["#d62728","#9467bd","#17becf"])

# #sns.stripplot(y = distances_to_plot, jitter=True, orient='h')
sns.swarmplot(data=df, x='Class', y='% of alternate structures\nwith this excision type',  orient='v', size=15, palette=palette)
# sns.swarmplot(y = fractional_distances_across_primary,  orient='h', size=10, color='b')
# sns.swarmplot(y = fractional_distances_inside_primary,  orient='h', size=10, color='g')

#n, bins, edges = plt.hist(distances_to_plot, bins= 100, density=True, alpha=0.5, label="Observed, mean=%2.2f" %(np.median(distances_to_plot)))

#n_alternate, bins_alternate, edges_alternate = plt.hist(distances_to_plot_alternate, bins= 20, density=True, alpha=0.5, label="Observed, mean=%2.2f" %(np.mean(distances_to_plot_alternate)))
plt.text(0, 0.3, 'mean= %2.2f\nmedian= %2.2f '%(np.mean(fractional_distances_inside_primary), np.median(fractional_distances_inside_primary)), fontsize=size-5, ha='center')
plt.text(1, 0.3, 'mean= %2.2f\nmedian= %2.2f '%(np.mean(fractional_distances_across_primary), np.median(fractional_distances_across_primary)), fontsize=size-5, ha='center')
plt.text(2, 0.3, 'mean= %2.2f\nmedian= %2.2f '%(np.mean(fractional_distances_shared_primary), np.median(fractional_distances_shared_primary)), fontsize=size-5, ha='center')

#plt.title("Regions near excision sites in primary structures are \n preferred in forming subsequent petite structures", fontsize=size)
plt.ylabel("Minimum fractional distance between\nprimary/alternate alignment edge", fontsize=size)
plt.xlabel("")
plt.yticks(fontsize=size, rotation=0)
#plt.text(1000, 0.0015, 'Observed: 66% of edges are within 580bp of an origin', fontsize=size)
#plt.text(0.2, 0.3, 'mean= %2.2f (1420bp)\nmedian= %2.2f (440bp) '%(np.mean(fractional_distances_to_plot), np.median(fractional_distances_to_plot)), fontsize=size)
#plt.ylabel("Frequency", fontsize=size)
#plt.xlabel("Source of junctions")
plt.xticks(fontsize=size, rotation=0)
#plt.yticks(fontsize=size, rotation=0)
#plt.xlim(0, 5000)
#plt.legend(loc='upper right', fontsize=size)
figure = plt.gcf()
figure.set_size_inches(14, 9)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

folder_name = os.sep + "ori_distance_dists" + os.sep
path = currDir + folder_name

if not os.path.exists(path):
    os.makedirs(path)

#save histogram figure
#plt.show()
plt.savefig(path + "primary_alternate_edge_distances_classes_std.png", bbox_inches="tight", dpi=300)

figure.clf()











