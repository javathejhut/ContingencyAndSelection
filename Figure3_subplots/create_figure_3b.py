import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random
import math
import vcf

origin_locations_forlines = [[4312, 4012], [32231, 32501], [54840, 54567], [56832, 56567], [82329, 82600] , [45227, 47927], [30220, 30594], [12780, 12510]]

#do our current data
selected_data = []
selected_simulated = []
selected_simulated_nonspanning = []

with open("pikl_pipeline_output_Normal_INREPEATS_data", 'r') as input_file:
    for line in input_file:
        selected_data.append(float(line.strip().split()[0]))

    selected_data = sorted(selected_data)

with open("pikl_pipeline_output_Normal_INREPEATS_simulated", 'r') as input_file:
    for line in input_file:
        selected_simulated.append(float(line.strip().split()[0]))

    selected_simulated = sorted(selected_simulated)

with open("pikl_pipeline_output_Normal_INREPEATS_simulated_nonspanning", 'r') as input_file:
    for line in input_file:
        selected_simulated_nonspanning.append(float(line.strip().split()[0]))

    selected_simulated_nonspanning = sorted(selected_simulated_nonspanning)

#file for blast located direct repeats

sample_name = "mtdna_self_deduplicated.out"

# blast dictionary
bD = {'QSEQID': 0, 'SSEQID': 1, 'PIDENT': 2, 'LENGTH': 3,
      'MISMATCH': 4, 'GAPOPEN': 5, 'QSTART': 6, 'QEND': 7,
      'SSTART': 8, 'SEND': 9, 'EVALUE': 10, 'BITSCORE': 11}

repeats = []
with open(sample_name, 'r') as blast_tab:
        for line in blast_tab:

            alignment_length = int(line.strip().split()[bD['LENGTH']])
            sstart = int(line.strip().split()[bD['SSTART']])
            send = int(line.strip().split()[bD['SEND']])
            qstart = int(line.strip().split()[bD['QSTART']])
            qend = int(line.strip().split()[bD['QEND']])
            bitscore = float(line.strip().split()[bD['BITSCORE']])
            pident = float(line.strip().split()[bD['PIDENT']])/100.0

            if send-sstart == qend-qstart:
                is_reversed = False
            else:
                is_reversed = True

            repeats.append([is_reversed, sstart, send, qstart, qend, alignment_length])

true_repeats = [ele for ele in repeats if ele[0]==False]
inverted_repeats = [ele for ele in repeats if ele[0]==True]

mean_true_repeat_locations = []
mean_inverted_repeat_locations = []
#compute mean locations of both true repeats and inverted repeats
for ele in true_repeats:
    mean_true_repeat_locations.append([int(np.mean([ele[1], ele[2]])), int(np.mean([ele[3], ele[4]])), ele[-1]])

for ele in inverted_repeats:
    mean_inverted_repeat_locations.append([int(np.mean([ele[1], ele[2]])), int(np.mean([ele[3], ele[4]])), ele[-1]])

delimeters = [23,22,21,20,19,18,17,16,15,14,13,12,11]

#
# for j, ele in enumerate(mean_true_repeat_locations):
#         if ele[-1] >= delim:
#             current_repeat_slice_true.append(ele)

distances_to_plot = []
for junction in mean_true_repeat_locations:
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

distances_to_plot = sorted(distances_to_plot)

#now make cumulatives:

#for selected
x =  np.arange(len(selected_data))
y = np.arange(len(selected_simulated))
z = np.arange(len(selected_simulated_nonspanning))
r = np.arange(len(distances_to_plot))

selected_data_cum = np.cumsum(np.ones_like(x)) / float(len(x))
selected_simulated_cum = np.cumsum(np.ones_like(y)) / float(len(y))
selected_simulated_nonspanning_cum = np.cumsum(np.ones_like(z)) / float(len(z))
repeats_cum = np.cumsum(np.ones_like(r)) / float(len(r))

print x
print selected_data
print selected_data_cum

for i in range(len(selected_data)):
    if selected_data[i] > 200:
        print selected_data[i]
	print selected_data_cum[i]
#now do plotting for all models
#plt.plot(lenient_data, lenient_data_cum, color = 'r', linestyle='solid', lw=2, label='lenient data')
#plt.plot(lenient_simulated, lenient_simulated_cum, color='r', linestyle='dotted', lw=2, label='lenient span simulated')

plt.plot(selected_data, selected_data_cum, color = 'blue', linestyle='solid', lw=3.5, label='Data')
plt.plot(selected_simulated, selected_simulated_cum, color='green', linestyle='dotted', lw=3.5, label='Simulated random alignments spanning origin')
plt.plot(selected_simulated_nonspanning, selected_simulated_nonspanning_cum, color='orange', linestyle='dashed', lw=3.5, label='Simulated random alignments')

plt.plot(distances_to_plot, repeats_cum, color='black', linestyle='dashed', lw=3.5, label='Perfect repeats on reference (>=11bp)')

# colours = ['C0', 'C1', 'C2', 'C3', 'C4']
# for i in range(3, 8):
#
#     plt.plot(chrI[i], cumulative_FP_chrI[i], color=colours[i-3], linestyle='solid', lw=2, label='%i-mer SSE model' % (i))
#     plt.plot(chrM[i], cumulative_FP_chrM[i], color=colours[i-3], linestyle='dotted', lw=2)

#plt.plot(chrI[model_number], cumulative_FP, linestyle='dotted',  lw=3, label='%i-mer SSE model'%(model_number))

size= 30
#plt.title("Cumulative distribution of Non-inverted junction\nedge distances to closest origins of replication", fontsize=size)
plt.xlabel("Displacement from closest origin edge (kbp)", fontsize=size)
plt.ylabel("Fraction of breakpoints", fontsize=size)
plt.yscale('log')
#plt.xscale('log')
plt.xlim(-2000,7000)
plt.xticks([-2000, -1000, 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000],[-2, -1, 0, 1, 2,3 ,4 ,5 ,6, 7], fontsize=size)
plt.ylim(0.001, 1.2)
plt.xticks(fontsize=size-0, rotation=0)
plt.yticks(fontsize=size-0, rotation=0)
plt.legend(loc='best', fontsize=size-9)

#plot section indicating inside repeats

line = np.arange(85779)

fill_negative = plt.fill_between(np.arange(-2000, 0), 0.001, 1.2 , color='gray', alpha= 0.3)

plt.text(-1850, 0.45, "Inside Origin", fontsize=size-10)

plt.text(4000, 0.45, "Outside Origin", fontsize=size-10)

figure = plt.gcf()
figure.set_size_inches(12, 10)

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir + os.sep

#path = currDir

#save histogram figure
#plt.show()
plt.savefig(path + "In_repeats_Normal.png", bbox_inches="tight", dpi=600)

#figure.clf()
#pvalue 10^-3 for 6mer, 10^-4 for 5mer for example

#now create cumulatives, with false positive SNP frequency for chrI control
# for key in chrI.keys():
#     x = np.arange(len(chrI[key]))
#     cumulative_FP = np.cumsum(np.ones_like(x))/float(len(x))
#     #print np.cumsum(chrI[key])[-10:-1]
#     plt.plot
#     print cumulative_FP[-10:]

# print chrI[3][-10:-1]
# print len(chrI[3])
# print len(chrI[4])


    #
    # plt.hist(counts[model_number], bins= 100, density=False, alpha=0.5, label='%i-mer SSE model'%(model_number))
    # #plt.xlim(0, 30000)
    # #plt.ylim(0, 2500)
    #
    # print model_number, min(counts[model_number])











    # vcf_reader=vcf.Reader(open(path_larges+ file, 'rb'))
    # if sample_name not in Grandes.keys():
    #     Grandes[sample_name] = []
    # for record in vcf_reader:
    #     if record.QUAL >= qual_threshold:
    #         if str(record.INFO['SVTYPE']) in genotyped_list:
    #             if str(record.INFO['SVTYPE']) == 'INV':
    #                 Grandes[sample_name].append([True, record.POS, record.INFO['END']])
    #
    #             else:
    #                 Grandes[sample_name].append([False, record.POS, record.INFO['END']])

# #store petites that pass filtering criteria
#
# for file in os.listdir(path_mediums):
#     sample_name = file.split('.')[0]
#
#     vcf_reader=vcf.Reader(open(path_mediums+ file, 'rb'))
#     if sample_name not in Petites.keys():
#         Petites[sample_name] = []
#     for record in vcf_reader:
#         if record.QUAL >= qual_threshold:
#             if str(record.INFO['SVTYPE']) in genotyped_list:
#                 if str(record.INFO['SVTYPE']) == 'INV':
#                     Petites[sample_name].append([True, record.POS, record.INFO['END']])
#
#                 else:
#                     Petites[sample_name].append([False, record.POS, record.INFO['END']])
#
# for file in os.listdir(path_smalls):
#     sample_name = file.split('.')[0]
#
#     vcf_reader=vcf.Reader(open(path_smalls+ file, 'rb'))
#     if sample_name not in Petites.keys():
#         Petites[sample_name] = []
#     for record in vcf_reader:
#         if record.QUAL >= qual_threshold:
#             if str(record.INFO['SVTYPE']) in genotyped_list:
#                 if str(record.INFO['SVTYPE']) == 'INV':
#                     Petites[sample_name].append([True, record.POS, record.INFO['END']])
#
#                 else:
#                     Petites[sample_name].append([False, record.POS, record.INFO['END']])
#
# #extracting normal junctions across all samples
# allsamples_normal_junctions = []
# allsamples_normal_junctions_E = []
#
# for sample in Grandes.keys():
#
#     if 'E' not in sample:
#         for junc in Grandes[sample]:
#             if not junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100 :
#                 allsamples_normal_junctions.append([junc[1], junc[2]])
#     else:
#         for junc in Grandes[sample]:
#             if not junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100:
#                 allsamples_normal_junctions_E.append([junc[1], junc[2]])
#
# for sample in Petites.keys():
#
#     if 'E' not in sample:
#         for junc in Petites[sample]:
#             if not junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100:
#                 allsamples_normal_junctions.append([junc[1], junc[2]])
#     else:
#         for junc in Petites[sample]:
#             if not junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100:
#                 allsamples_normal_junctions_E.append([junc[1], junc[2]])
#
# #extracting inverted junctions across all samples
# allsamples_inverted_junctions = []
# allsamples_inverted_junctions_E = []
#
# for sample in Grandes.keys():
#
#     if 'E' not in sample:
#         for junc in Grandes[sample]:
#             if junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100 :
#                 allsamples_inverted_junctions.append([junc[1], junc[2]])
#     else:
#         for junc in Grandes[sample]:
#             if junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100:
#                 allsamples_inverted_junctions_E.append([junc[1], junc[2]])
#
# for sample in Petites.keys():
#
#     if 'E' not in sample:
#         for junc in Petites[sample]:
#             if junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100:
#                 allsamples_inverted_junctions.append([junc[1], junc[2]])
#     else:
#         for junc in Petites[sample]:
#             if junc[0] and not abs(85779 - abs(junc[1]-junc[2])) <100:
#                 allsamples_inverted_junctions_E.append([junc[1], junc[2]])
#
#
# #print allsamples_normal_junctions
# #print allsamples_normal_junctions_E
#
# print len(allsamples_normal_junctions)
# print len(allsamples_normal_junctions_E)
# allsamples_normal_junctions.extend(allsamples_normal_junctions_E)
#
# #allsamples_inverted_junctions.extend(allsamples_inverted_junctions_E)
# print len(allsamples_normal_junctions)
#
# '''Extracting distances '''
#
# size=15
# distances_to_plot = []
# for junction in allsamples_normal_junctions:
#     distances = []
#     for ori in origin_locations_forlines:
#
#         #if contained in origin, consider this to be distance zero
#         if (int(junction[0]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))) or (
#                 int(junction[1]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))):
#
#             distances.append(0)
#         #if not contained, compute shortest distance
#         else:
#             dist1 = min(abs(int(junction[0]) - ori[0]), 85779 - abs(int(junction[0]) - ori[0]))
#             dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
#             dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
#             dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))
#
#             distances.append(min([dist1, dist2, dist3, dist4]))
#
#     distances_to_plot.append(min(distances))
#
# #do the same but for random distribution of junction locations
# '''Now do simulation based on these for specific repeat length delimeter'''
#
# #first compute length distribution of junctions
# normal_length_dist = []
# for junc in allsamples_normal_junctions:
#     normal_length_dist.append(int(abs(junc[1]-junc[0])))
#
# N_runs = 100
# min_distances_origin_simulated = []
# for i in range(N_runs):
#     simulated_junctions_per_run = []
#     for length in normal_length_dist:
#         start_loc = random.randint(0, 85779)
#         left_or_right = random.randint(0, 1)
#
#         if left_or_right:
#             end_loc = start_loc - length
#             if end_loc < 0:
#                 end_loc = end_loc + 85779
#             elif end_loc > 85779:
#                 end_loc = end_loc - 85779
#
#         else:
#             end_loc = start_loc + length
#             if end_loc < 0:
#                 end_loc = end_loc + 85779
#             elif end_loc > 85779:
#                 end_loc = end_loc - 85779
#
#         simulated_junctions_per_run.append([start_loc, end_loc])
#
#     for junction in simulated_junctions_per_run:
#
#         simulated_distances = []
#
#         for ori in origin_locations_forlines:
#             # if contained in origin, consider this to be distance zero
#             if (int(junction[0]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))) or (
#                     int(junction[1]) in range(min([ori[0], ori[1]]), max([ori[0], ori[1]]))):
#
#                 simulated_distances.append(0)
#             # if not contained, compute shortest distance
#             else:
#                 dist1 = min(abs(int(junction[0]) - ori[0]), 85779- abs(int(junction[0]) - ori[0]))
#                 dist2 = min(abs(int(junction[0]) - ori[1]), 85779 - abs(int(junction[0]) - ori[1]))
#                 dist3 = min(abs(int(junction[1]) - ori[0]), 85779 - abs(int(junction[1]) - ori[0]))
#                 dist4 = min(abs(int(junction[1]) - ori[1]), 85779 - abs(int(junction[1]) - ori[1]))
#
#                 simulated_distances.append(min([dist1, dist2, dist3, dist4]))
#
#         min_distances_origin_simulated.append(min(simulated_distances))
#
# #print distances_to_plot
# n, bins, edges = plt.hist(distances_to_plot, bins= 12, density=True, alpha=0.5, label="Observed, mean=%2.2f" %(np.mean(distances_to_plot)))
#
# # get the width of each bin
# bin_width = bins[1] - bins[0]
# bin_numbers_simulated = int(math.ceil((max(min_distances_origin_simulated) - min(min_distances_origin_simulated))/bin_width))
#
# print bin_width
# n_simulated, bins_simulated, edges_simulated = plt.hist(min_distances_origin_simulated, bins= bin_numbers_simulated, density=True, alpha=0.5, label="Random simulated, mean=%2.2f" %(np.mean(min_distances_origin_simulated)))
#
# # sum over number in each bin and mult by bin width, which can be factored out
# integral = bin_width * sum(n[0:1])
# integral_simulated = bin_width * sum(n_simulated[0:1])
#
# #compute percentage of distances below 580
# threshold = 600
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
# print perc_below, perc_below_simulated
#
# plt.title("Distances between alignment edges \nand closest origins of replication (Illumina)", fontsize=size)
# plt.xlabel("Distance between alignment edge and origin (bp)", fontsize=size)
# plt.text(2000, 0.0006, 'Observed: %2.2f percent of edges are within %ibp of an origin'%(perc_below*100, threshold), fontsize=size)
# plt.text(2000, 0.0005, 'Simulated: %2.2f percent of edges are within %ibp of an origin'%(perc_below_simulated*100, threshold), fontsize=size)
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
# plt.savefig(path + "normal_closest_ori_dist.png", bbox_inches="tight", DPI=300)
#
# figure.clf()
#
#
#









