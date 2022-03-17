import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import random
import math
import scipy
from scipy.optimize import curve_fit
import matplotlib as mpl

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)

mixed_samples = ["3I1", "4I1", "5I1","6I1", "10I1", "13I1", "17I1", "18I1", "20I1" , "23I1", "5I2",
  "8I1", "19I1"]

#read in suppressivity control values (10T3 x 1, which is grande x grande cross)
control_supp = []
control_supp_mean = 0
control_supp_std = 0
with open("control_supp_freq.txt", 'r') as supp_control_file:
    for line in supp_control_file:
        control_supp.append(float(line.strip().split()[1]))

control_supp_mean = np.mean(control_supp)
control_supp_std = np.std(control_supp)

#import from suppressivity file
suppressivity = {}
suppressivity_yerr = {}
suppressivity_name_to_progenitor = {}

with open("suppressivity.txt", 'r') as supp_file:
    for line in supp_file:
        name = line.strip().split('_')[0]
        supp = float(line.strip().split()[1])
        yerr = float(line.strip().split()[2])
        suppressivity[name] = supp

        suppressivity_yerr[name] = yerr
        suppressivity_name_to_progenitor[name] = line.strip().split()[0]

'''Gather all junctions/types in pickled output across all strains'''

# construct all dicts... for all types
alternate_repeat_lengths = {}
alternate_monomer_counts = {}

primary_repeat_lengths = {}
primary_repeat_lengths_double_mixed = {}

primary_monomer_counts = {}

# for petites
sorted_primary_edges = {}
sorted_alternate_edges = {}

for sample_pikl_file in glob.glob(os.path.normpath(os.path.join(os.getcwd(), os.pardir)) +'/pikl_pipeline_output_petites/*.pikl'):
    sample_name = sample_pikl_file.strip().split('/')[-1].split('_')[0]

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

        #print sample_name, sorted_primary_edges[sample_name]

        #compute primary repeat lengths from these
        repeat_length = 0
        for ele in sorted_primary_edges[sample_name]:
            repeat_length+= abs(int(ele[1])-int(ele[0]))

        primary_repeat_lengths[sample_name] = repeat_length

        if sample_name in mixed_samples:
            primary_repeat_lengths_double_mixed[sample_name] = 2* repeat_length
        else:
            primary_repeat_lengths_double_mixed[sample_name] = repeat_length

# ori1, 2, 3, .., 8
origin_locations_forlines = [[4312, 4012], [32231, 32501], [54840, 54567], [56832, 56567], [82329, 82600],
                             [45227, 47927], [30220, 30594], [12780, 12510]]
origin_locations_active = [[4312, 4012], [32231, 32501], [54840, 54567],
                           [82329, 82600]]  # ori1, 2, 3, 5 "active" according to bernardi


currDir = os.path.dirname(os.path.realpath(__file__))
folder_name = '/pikl_pipeline_output/'
path = currDir + folder_name

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

def return_count_density(query_list, reference_list):
    query_length = 0
    # perform pairwise query to reference alignment test
    ori_count = 0
    # print query_list
    for q_list in query_list:
        query_length += abs(min(q_list[0], q_list[1]) - max(q_list[0], q_list[1]))
        for r_list in reference_list:
            q_range = set(range(min(q_list[0], q_list[1]), max(q_list[0], q_list[1]) + 1))
            r_range = set(range(min(r_list[0], r_list[1]), max(r_list[0], r_list[1]) + 1))

            # symm_diff = (q_range ^ r_range) & q_range
            overlap = (q_range & r_range)
            if len(overlap)>0:
                ori_count +=1

    # print total_query_overlap_bp, query_length
    return float(ori_count) / query_length

# for alternates, return origin density
origin_fraction_primary = {}
origin_fraction_primary_active = {}
origin_count_density_primary = {}
origin_count_density_primary_active = {}

origin_fraction_wt_active = return_fractional_overlap([[0, 85779]], origin_locations_active)
origin_fraction_wt = return_fractional_overlap([[0, 85779]], origin_locations_forlines)

for key in sorted_primary_edges.keys():

    origin_fraction_primary[key] = [return_fractional_overlap(sorted_primary_edges[key], origin_locations_forlines)]
    origin_fraction_primary_active[key] = [
        return_fractional_overlap(sorted_primary_edges[key], origin_locations_active)]

    origin_count_density_primary[key] = [return_count_density(sorted_primary_edges[key],origin_locations_forlines)]
    origin_count_density_primary_active[key] = [return_count_density(sorted_primary_edges[key], origin_locations_active)]

#relevant parameters here will be length of primary repeat unit, and primary origin density
#plotting setup
size=32

inverted_primary = ["3I1", "4I1", "5I1", "6I1", "10I1", "13I1", "17I1", "18I1", "20I1", "23I1", "5I2", "8I1", "19I1"]

noninverted_primary = ["2I1", "10I2", "11I2", "12I1", "21I1", "22I1", "1I1", "7I1", "7I2", "8I2", "9I1", "9I2", "11I1", "12I2",
                       "13I2", "14I1", "15I1", "16I1", "20I2", "21I2", "22I2", "23I2", "24I1", "24I2", "6I2"]


print len(inverted_primary) + len(noninverted_primary)
families = ['_CP1', '_CP2', '_CP3', '_CP5', '_P2', '_P3', '_P4', '_P5', '_P7']
#family_colours = ['blue', 'green', 'red', 'cyan', 'maroon', 'orange', 'purple', 'yellow', 'magenta']
family_colours = ["#2CEE0E", "#F3715A", "#E3D200", "#5C2D91", "#ED1C24", "#72BF44", "#00FFFF", "#9D85BE", "#0369A3"]
supp_to_plot = []
supp_to_plot_nomixed = []

ori_frac_to_plot = []
ori_count_density_toplot = []
keys_supp_to_plot = []
keys_supp_to_plot_nomixed = []
repeat_lengths_to_plot = []
repeat_lengths_to_plot_double_mixed = []
for key in origin_fraction_primary.keys():
    #print primary_repeat_lengths[key]
    ori_frac_to_plot.append(origin_fraction_primary[key][0])
    ori_count_density_toplot.append(origin_count_density_primary[key][0])
    supp_to_plot.append(suppressivity[key])
    repeat_lengths_to_plot.append(primary_repeat_lengths[key])
    repeat_lengths_to_plot_double_mixed.append(primary_repeat_lengths_double_mixed[key])
    keys_supp_to_plot.append(key)

    if key=='12I1':
        print repeat_lengths_to_plot_double_mixed[-1], "testing length 12I1 here"

    if key not in mixed_samples:
        supp_to_plot_nomixed.append(suppressivity[key])
        keys_supp_to_plot_nomixed.append(key)

#most simple linear + exp theory
def theory_exp_linear(X, vGt, vPt):
    LG, LP, PG, PP = X

    num = LP * np.exp((vGt*PG - vPt*PP))
    denom = LG

    r = num/denom

    supp = 1.0/ (1.0 + r)

    return supp

def return_r_exp_linear(X, vGt, vPt):
    LG, LP, PG, PP = X

    num = LP * np.exp((vGt * PG - vPt * PP))
    denom = LG

    r = num / denom

    return r

#most exp theory (so no breakup)
def theory_exp(X, vGt, vPt):
    LG, LP, PG, PP = X

    num =  np.exp((vGt*PG - vPt*PP))
    denom = 1.0

    r = num/denom

    supp = 1.0/ (1.0 + r)

    return supp

def return_r_exp(X, vGt, vPt):
    LG, LP, PG, PP = X

    num =  np.exp((vGt * PG - vPt * PP))
    denom = 1.0

    r = num / denom

    return r

#linear theory
def theory_linear(X, vG_o_vP):
    LG, LP, PG, PP = X

    num =  LP*PG
    denom = LG*PP

    r = vG_o_vP*(num/denom)

    supp = 1.0/ (1.0 + r)

    return supp

def return_r_linear(X, vG_o_vP):
    LG, LP, PG, PP = X

    num = LP * PG
    denom = LG * PP

    r = vG_o_vP*(num/denom)

    return r

#linear theory no breakup (so no prefactors)
def theory_linear_nobreakup(X, vG_o_vP):
    LG, LP, PG, PP = X

    num =  PG
    denom = PP

    r = vG_o_vP*(num/denom)

    supp = 1.0/ (1.0 + r)

    return supp

def return_r_linear_nobreakup(X, vG_o_vP):
    LG, LP, PG, PP = X

    num = PG
    denom = PP

    r = vG_o_vP*(num/denom)

    return r


toplot_theorsupp = []

#arrays of values for specific subset of data we want to fit/plot
LG = [] #length of grande monomer
LP = [] #length of primary petite monomer
PG = [] #number density of origins in Grandes
PP = [] #number density of origins in primary monomer in petites

#populate arrays
grande_ori_count_density = 8.0/85779
grande_monomer_length = 85779.0

for key in origin_fraction_primary.keys():
    LG.append(grande_monomer_length)
    LP.append(primary_repeat_lengths_double_mixed[key])
    #LP.append(primary_repeat_lengths[key])
    PG.append(grande_ori_count_density)
    PP.append(origin_count_density_primary[key][0])

X = [LG, LP, PG, PP]
supp = np.asarray(supp_to_plot)

#now do fit for exp_linear with density
vGt = 1.0
vPt = 1.0
vG_o_vP = 1.0

models = ["exp_linear", "exp", "linear", "linear_nobreakup"]
models = ["exp_linear"]
for model in models:
    print model
    if model=="exp_linear":
        p0 = (vGt, vPt)

        popt, pcov = curve_fit(theory_exp_linear, X, supp, p0, bounds=(0.0, [10 ** 8, 10 ** 8]))
        theor_supp = theory_exp_linear(np.asarray(X), *popt)

        r_to_plot = return_r_exp_linear(np.asarray(X), *popt)
    elif model== "exp":
        p0 = (vGt, vPt)

        popt, pcov = curve_fit(theory_exp, X, supp, p0, bounds=(0.0, [10 ** 8, 10 ** 8]))
        theor_supp = theory_exp(np.asarray(X), *popt)

        r_to_plot = return_r_exp(np.asarray(X), *popt)

    elif model== "linear":
        p0 = (vG_o_vP)

        popt, pcov = curve_fit(theory_linear, X, supp, p0, bounds=(0.0, 10**8))
        theor_supp = theory_linear(np.asarray(X), *popt)

        r_to_plot = return_r_linear(np.asarray(X), *popt)

    elif model== "linear_nobreakup":
        p0 = (vG_o_vP)

        popt, pcov = curve_fit(theory_linear_nobreakup, X, supp, p0, bounds=(0.0, 10**8))
        theor_supp = theory_linear_nobreakup(np.asarray(X), *popt)

        r_to_plot = return_r_linear_nobreakup(np.asarray(X), *popt)

    # calculaing R squared value for plotting
    y = np.array(supp)
    yed = y.astype(np.float)

    yfit = np.array(theor_supp)
    yfited = yfit.astype(np.float)

    ss_res = np.sum((yed - yfited) ** 2)




    ss_tot = np.sum((yed - np.mean(yed)) ** 2)
    r_2 = 1 - (ss_res / ss_tot)

    print r_2, float(len(y)) * math.log(math.sqrt(2.0 * math.pi * ss_res / (len(y) - 1))) + (len(y) - 1) / 2.0
    '''plotting'''
    R = np.arange(0, 4.0, 0.1)
    supp_theor_func = lambda r: 1.0 / (1.0 + r)

    supp_theor_cont = np.array(map(supp_theor_func, R))

    variables = np.append(popt, r_2)

    if "exp" not in model:
        plt.plot(R, supp_theor_cont, 'k-', linewidth=3,
                 label=r"Equation (1) fit:""\n" r"$v_G/v_P= %4.0f$""\n" r"$R^2 = %2.2f$" % tuple(variables))


    else:

        plt.plot(R, supp_theor_cont, 'k-', linewidth=3,
                 label="Model")

    for i in range(len(supp_to_plot)):

        for c, fam in enumerate(families):
            if fam in suppressivity_name_to_progenitor[keys_supp_to_plot[i]]:
                colour=family_colours[c]
                #fam_name = suppressivity_name_to_progenitor[keys_supp_to_plot[i]].split('_')[1][:-2]
                fam_name = str(c+1)
                #print suppressivity_name_to_progenitor[keys_supp_to_plot[i]]

        plt.scatter(r_to_plot[i], supp_to_plot[i], color=colour, s=100, label='Petite family %s' %(fam_name))
        plt.errorbar(r_to_plot[i], supp_to_plot[i], suppressivity_yerr[keys_supp_to_plot[i]], color=colour, elinewidth=3)

    plt.xticks(fontsize=size, rotation=0)
    plt.yticks(fontsize=size, rotation=0)
    plt.ylim(0.1, 1.0)
    plt.xlabel(r"$N_G/N_P$", fontsize=size)
    plt.ylabel("Suppressivity", fontsize=size)
    # plt.show()

    figure = plt.gcf()
    figure.set_size_inches(12, 10)

    currDir = os.path.dirname(os.path.realpath(__file__))
    path = currDir + os.sep

    handles, labels = plt.gca().get_legend_handles_labels()

    by_label = OrderedDict(zip(labels, handles))

    #print by_label
    # plt.legend(by_label.values()[0:1] + by_label.values()[2:3] + by_label.values()[1:2] + by_label.values()[3:],
    #            by_label.keys()[0:1] + by_label.keys()[2:3] + by_label.keys()[1:2] + by_label.keys()[3:],
    #            fontsize=size - 10)

    plt.legend(by_label.values()[0:1] + by_label.values()[6:7] + by_label.values()[3:4] + by_label.values()[2:3] + by_label.values()[1:2] + by_label.values()[9:] + by_label.values()[7:8] +  by_label.values()[5:6] + by_label.values()[8:9] + by_label.values()[4:5],
               by_label.keys()[0:1] + by_label.keys()[6:7] + by_label.keys()[3:4] + by_label.keys()[2:3] + by_label.keys()[1:2] + by_label.keys()[9:] + by_label.keys()[7:8] +  by_label.keys()[5:6] + by_label.keys()[8:9] + by_label.keys()[4:5],
               fontsize=size - 10)

    # plt.legend(by_label.values()[0:1] + ,
    #            by_label.keys()[:] ,
    #            fontsize=size - 10)

    # plt.legend()
    #plt.savefig(path + "supp_func_r_%s.png"%(model), bbox_inches="tight", dpi=300, format='PNG')
    plt.savefig(path + "supp_func_r_%s.svg"%(model), bbox_inches="tight")

    plt.clf()



# #start off by plotting suppressivity as a function of R for non mixed samples
#
# R = np.arange(0, 3., 0.1)
# supp_theor_func = lambda r: 1.0/(1.0 +r)
#
# supp_theor_cont = np.array(map(supp_theor_func, R))
#
# print supp_theor_cont
# variables = np.append(popt_nomixed, r_2_nomixed)
# plt.plot(R, supp_theor_cont, 'k-', label=r"fit:""\n" r"$v_Gt^*= %4.3f$" "\n"r"$v_Pt^*=%5.3f$""\n" r"$R^2 = %2.2f$" % tuple(variables))
#
# for i in range(len(supp_to_plot_nomixed)):
#
#     for c, fam in enumerate(families):
#         if fam in suppressivity_name_to_progenitor[keys_supp_to_plot_nomixed[i]]:
#             colour=family_colours[c]
#             fam_name = suppressivity_name_to_progenitor[keys_supp_to_plot_nomixed[i]].split('_')[1][:-2]
#             #print suppressivity_name_to_progenitor[keys_supp_to_plot[i]]
#
#     plt.scatter(r_to_plot_nomixed[i], supp_to_plot_nomixed[i], color=colour, label='%s progeny' %(fam_name))
#     plt.errorbar(r_to_plot_nomixed[i], supp_to_plot_nomixed[i], suppressivity_yerr[keys_supp_to_plot_nomixed[i]], color=colour)
#
#     variables = np.append(popt_nomixed, r_2_nomixed)
#     variables_nodensity = np.append(popt_nodensity, r_2_nodensity)
#
#     # plt.scatter(r_to_plot_nomixed[i], theor_supp_nomixed[i], color='black', marker='x',label=r"fit:""\n" r"vG*t= %4.3f" "\n"r"vP*t=%5.3f""\n" r"$R^2 = %2.2f$" % tuple(variables))
#
#     # plt.scatter(r_to_plot_nodensity[i], theor_supp_nodensity[i], color='red', marker='x',
#     #             label=r"fit:""\n" r"vG*t= %4.3f" "\n"r"vP*t=%5.3f""\n" r"$R^2 = %2.2f$" % tuple(variables_nodensity))
#
#
#     # plt.scatter(ori_frac_to_plot[i], supp_to_plot[i], color=colour, label='Progeny of %s' % (fam_name))
#     # plt.errorbar(ori_frac_to_plot[i], supp_to_plot[i], suppressivity_yerr[keys_supp_to_plot[i]], color=colour)
#
# #show control supp
# #plt.plot(range(0, 85779), [control_supp_mean]*85779, 'gray', label='WT mating control')
# #plt.fill_between(range(0, 85779), [control_supp_mean + control_supp_std]*85779, [control_supp_mean - control_supp_std]*85779, alpha=0.5, color='gray')
# #plt.xlim(0, 85779)
# #plt.xticks([0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000], [0, 10, 20, 30, 40, 50, 60, 70, 80], fontsize=size, rotation=0)
# plt.xticks(fontsize=size, rotation=0)
# plt.yticks(fontsize=size, rotation=0)
# plt.ylim(0.1, 1.0)
# plt.xlabel(r"$R_G/R_P$", fontsize=size)
# plt.ylabel("Suppressivity", fontsize=size)
# #plt.show()
#
# figure = plt.gcf()
# figure.set_size_inches(12, 10)
#
# currDir = os.path.dirname(os.path.realpath(__file__))
# path = currDir + os.sep
#
# handles, labels = plt.gca().get_legend_handles_labels()
#
# by_label = OrderedDict(zip(labels, handles))
#
# plt.legend(by_label.values()[0:1] + by_label.values()[2:3] +by_label.values()[1:2] +by_label.values()[3:], by_label.keys()[0:1] + by_label.keys()[2:3] + by_label.keys()[1:2] + by_label.keys()[3:], fontsize=size-8)
# #plt.legend()
# plt.savefig(path + "supp_func_r_nomixed.png", bbox_inches="tight", dpi=300, format='PNG')
#
# plt.clf()

#start off by plotting suppressivity as a function of R for all samples
#
# R = np.arange(0, 4.0, 0.1)
# supp_theor_func = lambda r: 1.0/(1.0 +r)
#
# supp_theor_cont = np.array(map(supp_theor_func, R))
# print popt
# variables = np.append(popt, r_2)
# plt.plot(R, supp_theor_cont, 'k-', label=r"fit:""\n" r"$v_Gt^*= %4.3f$" "\n"r"$v_Pt^*=%5.3f$""\n" r"$R^2 = %2.2f$" % tuple(variables))
#
# for i in range(len(supp_to_plot)):
#
#     for c, fam in enumerate(families):
#         if fam in suppressivity_name_to_progenitor[keys_supp_to_plot[i]]:
#             colour=family_colours[c]
#             fam_name = suppressivity_name_to_progenitor[keys_supp_to_plot[i]].split('_')[1][:-2]
#             #print suppressivity_name_to_progenitor[keys_supp_to_plot[i]]
#
#     plt.scatter(r_to_plot[i], supp_to_plot[i], color=colour, label='%s progeny' %(fam_name))
#     plt.errorbar(r_to_plot[i], supp_to_plot[i], suppressivity_yerr[keys_supp_to_plot[i]], color=colour)
#
#     # plt.scatter(r_to_plot_nomixed[i], theor_supp_nomixed[i], color='black', marker='x',label=r"fit:""\n" r"vG*t= %4.3f" "\n"r"vP*t=%5.3f""\n" r"$R^2 = %2.2f$" % tuple(variables))
#
#     # plt.scatter(r_to_plot_nodensity[i], theor_supp_nodensity[i], color='red', marker='x',
#     #             label=r"fit:""\n" r"vG*t= %4.3f" "\n"r"vP*t=%5.3f""\n" r"$R^2 = %2.2f$" % tuple(variables_nodensity))
#
#
#     # plt.scatter(ori_frac_to_plot[i], supp_to_plot[i], color=colour, label='Progeny of %s' % (fam_name))
#     # plt.errorbar(ori_frac_to_plot[i], supp_to_plot[i], suppressivity_yerr[keys_supp_to_plot[i]], color=colour)
#
# #show control supp
# #plt.plot(range(0, 85779), [control_supp_mean]*85779, 'gray', label='WT mating control')
# #plt.fill_between(range(0, 85779), [control_supp_mean + control_supp_std]*85779, [control_supp_mean - control_supp_std]*85779, alpha=0.5, color='gray')
# #plt.xlim(0, 85779)
# #plt.xticks([0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000], [0, 10, 20, 30, 40, 50, 60, 70, 80], fontsize=size, rotation=0)
# plt.xticks(fontsize=size, rotation=0)
# plt.yticks(fontsize=size, rotation=0)
# plt.ylim(0.1, 1.0)
# plt.xlabel(r"$R_G/R_P = \frac{1}{L_G}e^{V_G \rho_G t^*}/\frac{1}{L_P}e^{V_P \rho_P t^*}$", fontsize=size)
# plt.ylabel("Suppressivity", fontsize=size)
# #plt.show()
#
# figure = plt.gcf()
# figure.set_size_inches(12, 10)
#
# currDir = os.path.dirname(os.path.realpath(__file__))
# path = currDir + os.sep
#
# handles, labels = plt.gca().get_legend_handles_labels()
#
# by_label = OrderedDict(zip(labels, handles))
#
# plt.legend(by_label.values()[0:1] + by_label.values()[2:3] +by_label.values()[1:2] +by_label.values()[3:], by_label.keys()[0:1] + by_label.keys()[2:3] + by_label.keys()[1:2] + by_label.keys()[3:], fontsize=size-10)
# #plt.legend()
# plt.savefig(path + "supp_func_r.png", bbox_inches="tight", dpi=300, format='PNG')
#
# plt.clf()
#
# #for ori count density
# for i in range(len(supp_to_plot)):
#
#     for c, fam in enumerate(families):
#         if fam in suppressivity_name_to_progenitor[keys_supp_to_plot[i]]:
#             colour=family_colours[c]
#             fam_name = suppressivity_name_to_progenitor[keys_supp_to_plot[i]].split('_')[1][:-2]
#
#     plt.scatter(ori_count_density_toplot[i], supp_to_plot[i], color=colour, label='%s progeny' %(fam_name))
#     plt.errorbar(ori_count_density_toplot[i], supp_to_plot[i], suppressivity_yerr[keys_supp_to_plot[i]], color=colour)
#
#     #plt.scatter(ori_frac_to_plot[i], supp_to_plot[i], color=colour, label='Progeny of %s' % (fam_name))
#     #plt.errorbar(ori_frac_to_plot[i], supp_to_plot[i], suppressivity_yerr[keys_supp_to_plot[i]], color=colour)
#
# #also include a line indicating WT ori count density
# plt.plot([8./85779]*11, np.arange(0, 1.1, 0.1), 'g--', linewidth=2)
#
# plt.xticks(fontsize=size, rotation=0)
# plt.yticks(fontsize=size, rotation=0)
# #plt.ylim(0, 1.0)
# plt.xlabel("Origin count density in monomer", fontsize=size)
# plt.ylabel("Suppressivity", fontsize=size)
# #plt.show()
# plt.xlim(0, 0.0012)
# figure = plt.gcf()
# figure.set_size_inches(12, 10)
#
# currDir = os.path.dirname(os.path.realpath(__file__))
# path = currDir + os.sep
#
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = OrderedDict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys(), fontsize=size-5)
#
# plt.savefig(path + "suppressivity_vs_origin_count_density.png", bbox_inches="tight", dpi=300, format='PNG')
#
# plt.clf()

















