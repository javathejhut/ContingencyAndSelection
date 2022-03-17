import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random
import math

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

        keys_in_repeats_normal = []
        keys_in_repeats_inverted = []
        for key in junction_counts_dict.keys():
            if junction_type_dict[key] == 'False':
                if key in primary_aln_key_count[0] and key not in keys_in_repeats_normal:
                    keys_in_repeats_normal.append(key)

                else:
                    for ele in alternate_aln_key_count:
                        if key in ele[0] and key not in keys_in_repeats_normal:
                            keys_in_repeats_normal.append(key)
            else:
                if key in primary_aln_key_count[0] and key not in keys_in_repeats_inverted:
                    keys_in_repeats_inverted.append(key)

                else:
                    for ele in alternate_aln_key_count:
                        if key in ele[0] and key not in keys_in_repeats_inverted:
                            keys_in_repeats_inverted.append(key)

        # now extract repeats that have these
        # junctions_in_repeats[sample_name] = [
        #     [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in
        #     keys_in_repeats]

        allsamples_normal_junctions.extend([
            [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in
            keys_in_repeats_normal])

        allsamples_inverted_junctions.extend([
            [junction_ref_boundaries_dict[key][0][2], junction_ref_boundaries_dict[key][1][2]] for key in
            keys_in_repeats_inverted])

#
# allsamples_normal_junctions.remove([1, 85777])
# allsamples_normal_junctions.remove([1, 85778])
# allsamples_normal_junctions.remove([1, 85776])
# allsamples_normal_junctions.remove([20, 85770])
# allsamples_normal_junctions.remove([48, 85750])

'''Extracting distances '''

if junction_type=='Normal':
    allsamples_normal_junctions = allsamples_normal_junctions
else:
    allsamples_normal_junctions = allsamples_inverted_junctions

size=15
distances_to_plot = []
for junction in allsamples_normal_junctions:
    distances = []
    for ori in origin_locations_forlines:

        # if contained in origin, consider this to be distance negative
        if (int(junction[0]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))) or (
                int(junction[1]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))):

            dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
            dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
            dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
            dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))

            distances.append(min([-abs(dist1), -abs(dist2), -abs(dist3), -abs(dist4)], key=abs))
        # if not contained, compute shortest distance
        else:
            dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
            dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
            dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
            dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))

            distances.append(min([dist1, dist2, dist3, dist4]))

    distances_to_plot.append(min(distances, key=abs))

'''Now do simulation INCLUDING ORIGINS based on these for specific repeat length delimeter'''

#first compute length distribution of junctions
normal_length_dist = []
for junc in allsamples_normal_junctions:
    normal_length_dist.append(int(abs(junc[1]-junc[0])))

N_runs = 1000
min_distances_spanning_origins_simulated = []

# print normal_length_dist
# print np.mean(normal_length_dist)
# print np.median(normal_length_dist)
for i in range(N_runs):

    simulated_junctions_per_run = []
    for length in normal_length_dist:

        #first pick origin at random
        origin_index = random.randint(0, 7)

        origin_loc_selected = [min(origin_locations_forlines[origin_index]), max(origin_locations_forlines[origin_index])]

        if origin_loc_selected[0] < length:
            start_range = 0
            end_range = -(origin_loc_selected[0] - length) + origin_loc_selected[1]

            start_loc_add = random.randint(start_range, end_range)

            start_loc = 85779 + (origin_loc_selected[0] - length) +start_loc_add

            if start_loc > 85779:
                start_loc = start_loc - 85779

            end_loc = start_loc + length

        else:

            start_loc = random.randint(origin_loc_selected[0] - length, origin_loc_selected[1])

        end_loc = start_loc + length

        if end_loc > 85779:
            end_loc = end_loc - 85779

        #print start_loc, end_loc, length

        simulated_junctions_per_run.append([start_loc, end_loc])

    for junction in simulated_junctions_per_run:

        simulated_distances = []

        for ori in origin_locations_forlines:

            # if contained in origin, consider this to be distance negative
            if (int(junction[0]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))) or (
                    int(junction[1]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))):

                dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
                dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
                dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
                dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))

                simulated_distances.append(min([-abs(dist1), -abs(dist2), -abs(dist3), -abs(dist4)], key=abs))
            # if not contained, compute shortest distance
            else:
                dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
                dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
                dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
                dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))

                simulated_distances.append(min([dist1, dist2, dist3, dist4]))

        min_distances_spanning_origins_simulated.append(min(simulated_distances, key=abs))


#do the same but for random distribution of junction locations
'''Now do simulation based on these for specific repeat length delimeter'''

#first compute length distribution of junctions
normal_length_dist = []
for junc in allsamples_normal_junctions:
    normal_length_dist.append(int(abs(junc[1]-junc[0])))

N_runs = 1000
min_distances_origin_simulated = []
for i in range(N_runs):
    simulated_junctions_per_run = []
    for length in normal_length_dist:
        start_loc = random.randint(0, 85779)
        left_or_right = random.randint(0, 1)

        if left_or_right:
            end_loc = start_loc - length
            if end_loc < 0:
                end_loc = end_loc + 85779
            elif end_loc > 85779:
                end_loc = end_loc - 85779

        else:
            end_loc = start_loc + length
            if end_loc < 0:
                end_loc = end_loc + 85779
            elif end_loc > 85779:
                end_loc = end_loc - 85779

        simulated_junctions_per_run.append([start_loc, end_loc])

    for junction in simulated_junctions_per_run:

        simulated_distances = []

        for ori in origin_locations_forlines:
            # if contained in origin, consider this to be distance negative
            if (int(junction[0]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))) or (
                    int(junction[1]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))):

                dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
                dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
                dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
                dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))

                simulated_distances.append(min([-abs(dist1), -abs(dist2), -abs(dist3), -abs(dist4)], key=abs))
            # if not contained, compute shortest distance
            else:
                dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
                dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
                dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
                dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))

                simulated_distances.append(min([dist1, dist2, dist3, dist4]))

        min_distances_origin_simulated.append(min(simulated_distances, key=abs))

#print distances_to_plot
fixed_bins = np.arange(-2000, 13000, 500)
n, bins, edges = plt.hist(distances_to_plot, bins= fixed_bins, density=True, alpha=0.5, label="Observed, mean=%2.2f" %(np.mean(distances_to_plot)))

# get the width of each bin
bin_width = bins[1] - bins[0]
bin_numbers_simulated = int(math.ceil((max(min_distances_origin_simulated) - min(min_distances_origin_simulated))/bin_width))

n_simulated, bins_simulated, edges_simulated = plt.hist(min_distances_origin_simulated, bins= fixed_bins, density=True, alpha=0.5, label="Random simulated, mean=%2.2f" %(np.mean(min_distances_origin_simulated)))
n_simulated_span, bins_simulated_span, edges_simulated_span = plt.hist(min_distances_spanning_origins_simulated, bins=fixed_bins, density=True, alpha=0.5, label="Random simulated span origins, mean=%2.2f" %(np.mean(min_distances_spanning_origins_simulated)))

print len(distances_to_plot)
print len(min_distances_spanning_origins_simulated)

with open(currDir + '/ori_distance_dists/' + folder_name[1:-1] + "_" + str(junction_type) +"_INREPEATS_simulated", 'w') as output_simulated:
    for ele in min_distances_spanning_origins_simulated:
        output_simulated.write(str(ele) + "\n")

with open(currDir + '/ori_distance_dists/' + folder_name[1:-1] + "_" +str(junction_type) +"_INREPEATS_data", 'w') as output_data:
    for ele in distances_to_plot:
        output_data.write(str(ele) + "\n")

with open(currDir + '/ori_distance_dists/' + folder_name[1:-1] + "_" +str(junction_type) + "_INREPEATS_simulated_nonspanning", 'w') as output_data:
    for ele in min_distances_origin_simulated:
        output_data.write(str(ele) + "\n")

# # sum over number in each bin and mult by bin width, which can be factored out
# integral = bin_width * sum(n[0:4])
# integral_simulated = bin_width * sum(n_simulated[0:4])
#
# #compute percentage of distances below 580
# threshold = 0
# total_count = 0.0
# below_thresh_count = 0.0
# for dist in distances_to_plot:
#     total_count+=1.0
#     if dist <=threshold:
#         below_thresh_count+=1.0
#
# perc_below = below_thresh_count/total_count
#
# total_count = 0.0
# below_thresh_count = 0.0
# for dist in min_distances_origin_simulated:
#     total_count+=1.0
#     if dist <=threshold:
#         below_thresh_count+=1.0
#
# perc_below_simulated = below_thresh_count/total_count
#
# total_count = 0.0
# below_thresh_count = 0.0
# for dist in min_distances_spanning_origins_simulated:
#     total_count+=1.0
#     if dist <=threshold:
#         below_thresh_count+=1.0
#
# perc_below_simulated_span = below_thresh_count/total_count
#
# print perc_below, perc_below_simulated, perc_below_simulated_span
#
# print bins[0:4]
# print integral, integral_simulated
#
#
# plt.title("Distances between alignment edges \nand closest edges of origins of replication (Nanopore)", fontsize=size)
# plt.xlabel("Distance between alignment edge and origin (bp)", fontsize=size)
# plt.text(2000, 0.00045, 'Observed: %2.2f percent of edges are within an origin \n(-ve distance)'%(perc_below*100), fontsize=size)
# plt.text(2000, 0.0004, 'Simulated: %2.2f percent of edges are within an origin \n(-ve distance)'%(perc_below_simulated*100), fontsize=size)
# plt.text(2000, 0.00035, 'Simulated spanning origins: %2.2f percent of edges are \n within an origin (-ve distance)'%(perc_below_simulated_span*100), fontsize=size)
# plt.ylabel("Frequency", fontsize=size)
# plt.xticks(fontsize=size, rotation=0)
# plt.yticks(fontsize=size, rotation=0)
# plt.legend(loc='upper right', fontsize=size)
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
# plt.savefig(path + "normal_closest_ori_dist_incl_originspan.png", bbox_inches="tight", DPI=300)
#
# figure.clf()











