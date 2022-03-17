import pysam
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import statistics
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
import pickle

######GLOBAL PARAMS#######
del_ins_threshold = 30 #read space vs ref space mismatch required to call junction
min_alignment_length = 300 #mininmum number of bp for alignment to have (removed before junction calling)
min_read_threshold = 3 #minimum number of independent reads that must support one cluster
cluster_merging_eps = 1000 #maximum distance between points to cluster in junction space (note invdup is automatic eps)
binomial_merging_frac = 0.34 #fractional difference from expected mean of pattern below which junction transitions are merged

#print invdup_loc_prob
def return_prob_invdup_location_artifact(frac_location):
    for i in range(len(invdup_readloc_bins) -1):
        start_interval = invdup_readloc_bins[i]
        end_interval = invdup_readloc_bins[i+1]

        if (frac_location <= end_interval) and (frac_location >= start_interval):
            return invdup_loc_prob[i]

def complimentary_strand(a_string):
    to_return = ""

    for c in a_string:
        to_return += complimentary[c]

    return to_return[::-1]

def is_fully_overlapped(ref_start, ref_end, query_start, query_end):
    query_range = [int(query_start), int(query_end)]
    if (int(ref_start) in range(min(query_range), max(query_range))) \
            and (int(ref_end) in range(min(query_range), max(query_range))):
        return True
    else:
        return False

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

def return_junctions_and_locations(subalignments, deletion_threshold):
    #recall, subaligments take form [read start, read end, ref start, ref end, reversed=True]

    junctions = []
    junction_read_locations = []
    endpoints = [] #start/end reference mapping position with orientations
    deletion_insertion_threshold = deletion_threshold #minimum deviation in read vs reference coordinates to call a SV

    if len(subalignments)>0:
        endpoints = [[subalignments[0][2],subalignments[0][4]],[subalignments[-1][3], subalignments[-1][4]]]

    for i in range(len(subalignments)-1):
        curr_read_end = subalignments[i][1]
        next_read_start = subalignments[i+1][0]

        curr_ref_end = subalignments[i][3]
        next_ref_start = subalignments[i+1][2]

        curr_orientation = subalignments[i][4]
        next_orientation = subalignments[i+1][4]

        name = subalignments[i][6]

        #options are insertions/deletions
        if curr_orientation == next_orientation:

            #deletions
            if (abs(next_ref_start - curr_ref_end) - abs(next_read_start - curr_read_end)) > deletion_insertion_threshold:
                junctions.append([max(curr_ref_end, next_ref_start),min(curr_ref_end, next_ref_start) , False, curr_orientation, next_orientation, curr_ref_end, next_ref_start, name])
                junction_read_locations.append(statistics.mean([curr_read_end,next_read_start]))

            #insertions
            elif abs(next_read_start - curr_read_end) - abs(next_ref_start - curr_ref_end) > deletion_insertion_threshold:
                junctions.append([max(curr_ref_end, next_ref_start),min(curr_ref_end, next_ref_start) , False, curr_orientation, next_orientation, curr_ref_end, next_ref_start,name])
                junction_read_locations.append(statistics.mean([curr_read_end, next_read_start]))

        elif curr_orientation != next_orientation:

            junctions.append([min(curr_ref_end, next_ref_start), max(curr_ref_end, next_ref_start), True, curr_orientation, next_orientation,curr_ref_end, next_ref_start, name])
            junction_read_locations.append(statistics.mean([curr_read_end, next_read_start]))

            # #if we have an overlap, return average reference location at this position
            # if (next_read_start - curr_read_end) < 0:
            #     avg_ref_location = (curr_ref_end + next_ref_start)/2.0
            #     junctions.append([avg_ref_location, avg_ref_location, True])
            # else:
            #     junctions.append([min(curr_ref_end, next_ref_start), max(curr_ref_end, next_ref_start), True])

    return junctions, junction_read_locations, endpoints

'''For plausible inverted duplications'''

def compute_optimal_eps(class_type):
    # first trying to determine "optimal" epsilon algorithmically
    juncs_in_coords = np.array([[ele[0], ele[1]] for ele in class_type])
    neigh = NearestNeighbors(n_neighbors=2)
    nbrs = neigh.fit(juncs_in_coords)
    distances, indices = nbrs.kneighbors(juncs_in_coords)

    distances = np.sort(distances, axis=0)
    distances = distances[:, 1]

    pt1 = (len(distances) - 1, distances[-1])
    pt2 = (0, distances[0])

    diff_distances = []

    for i in range(0, len(distances)):
        pt = (i, distances[i])
        v21 = np.subtract(pt2, pt1)
        v10 = np.subtract(pt1, pt)
        distance = np.linalg.norm(np.cross(v21, v10)) / np.linalg.norm(v21)
        diff_distances.append(distance)
    ind = np.argmax(diff_distances)
    return distances[ind]

def compute_DBSCAN(our_eps, minsamples, class_type):
    juncs_in_coords = np.array([[ele[0], ele[1]] for ele in class_type])
    db = DBSCAN(eps=our_eps, min_samples=minsamples).fit(juncs_in_coords)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    return labels, core_samples_mask

def count_cluster_sizes(labels, core_samples_mask, class_type):
    class_counts = []
    for k in set(labels):

        class_member_mask = (labels == k)

        #xy = class_type[class_member_mask & core_samples_mask]
        xy_total_class_counting = np.array(class_type)[class_member_mask]

        #if not part of the noise cluster
        if k!=-1:
            class_counts.append(len(xy_total_class_counting))
    return class_counts

def plot_class_type(labels, class_type, core_samples_mask, scalar_map, string_type):
    juncs_in_coords = np.array([[int(ele[0]), int(ele[1])] for ele in class_type])

    unique_labels = set(labels)
    edge_colours = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

    for k, edge_col in zip(unique_labels, edge_colours):

        class_member_mask = (labels == k)

        xy = juncs_in_coords[class_member_mask & core_samples_mask]
        count = len(np.array(class_type)[class_member_mask])

        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
        else:
            col = scalar_map.to_rgba(count)
        if string_type=="invdup" or string_type=="inverted":
            marker='o'
            print "here"
        elif string_type == "normal":
            marker='s'
        if len(xy[:, 0]) > 0 and len(xy[:, 1]) > 0:
            # mean_x = np.mean(xy[:, 0])
            # mean_y = np.mean(xy[:, 1])
            # if not math.isnan(mean_x) and not math.isnan(mean_y):
            #     class_centroids.append([mean_x, mean_y])

            plt.plot(xy[:, 0], xy[:, 1], marker, markerfacecolor=col,
                     markeredgecolor=tuple(edge_col), markersize=10, alpha=0.4,  markeredgewidth=2)

        xy = juncs_in_coords[class_member_mask & ~core_samples_mask]

        if len(xy[:, 0]) > 0 and len(xy[:, 1]) > 0:
            plt.plot(xy[:, 0], xy[:, 1], marker, markerfacecolor=col,
                     markeredgecolor=tuple(edge_col), markersize=6, alpha=0.4)

def plot_artifact_noise(class_type):
    juncs_in_coords = np.array([[int(ele[0]), int(ele[1])] for ele in class_type])
    col = [0, 0, 0, 1]

    marker = 'o'
    xy = juncs_in_coords
    if len(xy[:, 0]) > 0 and len(xy[:, 1]) > 0:
        plt.plot(xy[:, 0], xy[:, 1], marker, markerfacecolor=col,
                 markeredgecolor='k', markersize=6, alpha=0.4)

def filter_cluster_read_support(labels, class_type, core_samples_mask, min_read_support):
    #juncs_in_coords = np.array([[int(ele[0]), int(ele[1])] for ele in class_type])
    juncs_array = np.array(class_type)
    unique_labels = set(labels)

    labels_to_remove = []
    readnames_per_cluster = []
    points_in_cluster = []
    unique_readnames_count = 0
    for k in unique_labels:

        if k!=-1:

            class_member_mask = (labels == k)
            points_in_cluster = juncs_array[class_member_mask]
            #record readnames contributing to cluster, then take set to count
            readnames_per_cluster = [junc[-1] for junc in points_in_cluster]
            unique_readnames_count = len(set(readnames_per_cluster))

            if unique_readnames_count < min_read_support:
                labels_to_remove.append(k)

    for removed in labels_to_remove:
        for n, ele in enumerate(labels):
            if ele==removed:
                labels[n] = -1
    return labels, core_samples_mask

def rotate(some_list):
    to_return = []
    length = len(some_list)
    for i in range(length):
        to_return.append(some_list[-i:] + some_list[:-i])
        to_return.append([-1* ele for ele in to_return[-1]])
    return to_return

def permute_only(some_list):
    to_return = []
    length = len(some_list)
    for i in range(length):
        to_return.append(some_list[-i:] + some_list[:-i])
        #to_return.append([-1* ele for ele in to_return[-1]])
    return to_return


def lcp(s,t):
    n = min(len(s), len(t))
    for i in range(0,n):
        if s[i]!=t[i]:
            return s[0:i]
    return s[0:n]

def lrs(read_list):
    lrs = []
    n = len(read_list)
    for i in range(0,n):
        for j in range(i+1, n):
            x = lcp(read_list[i:n], read_list[j:n])

            if len(x) > len(lrs):
                    lrs = x

    return lrs

def allrs(read_list):
    rs_dict = {}
    n = len(read_list)
    for i in range(0,n):
        for j in range(i+1, n):
            x = lcp(read_list[i:n], read_list[j:n])

            x_tup = tuple(x)
            if len(x_tup)!=0 and x_tup not in rs_dict.keys():
                rs_dict[x_tup] = 0

    for key in rs_dict.keys():
        length = len(key)
        #traverse original list in sliding window of current repeated string

        for j in range(0, len(read_list) - length + 1):
            window = tuple(read_list[j:j+length])
            if window == key:
                rs_dict[key]+=1

    return rs_dict

def allrs_indicies(read_list):
    rs_dict = {}
    rs_loc_dict = {}
    n = len(read_list)
    for i in range(0,n):
        for j in range(i+1, n):
            x = lcp(read_list[i:n], read_list[j:n])

            x_tup = tuple(x)
            if len(x_tup)!=0 and x_tup not in rs_dict.keys():
                rs_dict[x_tup] = 0
                rs_loc_dict[x_tup] = []

    for key in rs_dict.keys():
        length = len(key)
        #traverse original list in sliding window of current repeated string

        for j in range(0, len(read_list) - length + 1):
            rs_start = j
            rs_end = j + length

            window = tuple(read_list[j:j+length])
            if window == key:
                rs_dict[key]+=1
                rs_loc_dict[key].append([rs_start, rs_end])

    return rs_dict, rs_loc_dict

def return_l_mf_repeats(repeat_dict):
    repeat_lengths = []
    max_len_repeats = []
    max_len_repeat_freq = []
    repeats_to_return = []

    for key in repeat_dict:
        repeat_lengths.append(len(key))

    if len(repeat_lengths) > 0:
        max_length = max(repeat_lengths)

        #record repeats that correspond to "maximum" length
        for key in repeat_dict:
            if len(key) == max_length:
                max_len_repeats.append(key)

        #if we have multiple maximum length repeats,
        #return repeats with maximum frequency
        for repeat in max_len_repeats:
            max_len_repeat_freq.append(repeat_dict[repeat])

        max_len_max_freq = max(max_len_repeat_freq)

        for repeat in max_len_repeats:
            if repeat_dict[repeat] == max_len_max_freq:
                repeats_to_return.append(repeat)

        return repeats_to_return
    else:
        return []

def return_l_mf_repeats_indicies(repeat_dict, repeat_loc_dict):
    repeat_lengths = []
    max_len_repeats = []
    max_len_repeat_freq = []
    repeats_to_return = []
    repeat_loc_to_return = []

    for key in repeat_dict:
        repeat_lengths.append(len(key))

    if len(repeat_lengths) > 0:
        max_length = max(repeat_lengths)

        #record repeats that correspond to "maximum" length
        for key in repeat_dict:
            if len(key) == max_length:
                max_len_repeats.append(key)

        #if we have multiple maximum length repeats,
        #return repeats with maximum frequency
        for repeat in max_len_repeats:
            max_len_repeat_freq.append(repeat_dict[repeat])

        max_len_max_freq = max(max_len_repeat_freq)

        for repeat in max_len_repeats:
            if repeat_dict[repeat] == max_len_max_freq:
                repeats_to_return.append(repeat)
                repeat_loc_to_return.append(repeat_loc_dict[repeat]) #all locations for max repeat

        return repeats_to_return, repeat_loc_to_return
    else:
        return [], []

def uniform(read_list):
    occupants = []
    for ele in read_list:
        if ele not in occupants:
            occupants.append(ele)

    if len(occupants) == 1:
        return True
    else:
        return False

def update_vote_dict(test, the_dict):
    letters_with_positive_votes = []
    letters_with_negative_votes = []

    # for reads with no artifacts, return LRSes
    if '-' not in test:
        test_abs = [abs(ele) if isinstance(ele, int) else ele for ele in test]
        letters_present = list(set(test_abs))

        #update votes
        positive_votes =  return_l_mf_repeats(allrs(test_abs))

        if len(positive_votes) >0:

            for tup in positive_votes:
                for num in tup:
                    letters_with_positive_votes.append(num)

            letters_with_positive_votes = list(set(letters_with_positive_votes))
            letters_with_negative_votes = [lett for lett in letters_present if lett not in letters_with_positive_votes]

        else:
            letters_with_positive_votes = []
            letters_with_negative_votes = []

        #update votes

        for letter in letters_with_negative_votes:
            if letter != '*' and letter != 0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][1]+=1
                else:
                    the_dict[letter][1] += 1

        for letter in letters_with_positive_votes:
            if letter != '*' and letter != 0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][0]+=1
                else:
                    the_dict[letter][0] += 1

    #for reads with artifacts, return any abs() repeated characters across aligned sections
    else:

        indicies_delim = [i for i, n in enumerate(test) if n=='-']
        alignment_blocks = [] #list of set of chars between asterisks
        test_abs = [abs(ele) if isinstance(ele, int) else ele for ele in test]
        letters_present = list(set(test_abs))

        if len(indicies_delim)==1:
            #handle start:
            alignment_blocks.append(list(test_abs[0: indicies_delim[0]]))

            #handle end:
            alignment_blocks.append(list(test_abs[indicies_delim[0] + 1:]))

            #do intersection of both halves////////

        elif len(indicies_delim)>1:
            # handle start sections:
            alignment_blocks.append(list(test_abs[0: indicies_delim[0]]))

            # handle middle sections
            mid_ind = indicies_delim[0:len(indicies_delim)]
            for i in range(len(mid_ind)-1):
                alignment_blocks.append(list(test_abs[mid_ind[i] + 1: mid_ind[i+1]]))

            # handle end sections
            alignment_blocks.append(list(test_abs[indicies_delim[-1] + 1:]))

            #do pairwise intersection of all sections and record, then deduplicate

        #print alignment_blocks
        intersections = []
        for i, junc_list in enumerate(alignment_blocks):
            for j, junc_list2 in enumerate(alignment_blocks):
                if i!=j:
                    intersections.append(list(set(junc_list) & set(junc_list2)))

        #print intersections
        positive_votes = [list(item) for item in set(tuple(row) for row in intersections)]
        letters_with_positive_votes = []
        for some_list in positive_votes:
            for ele in some_list:
                if ele not in letters_with_positive_votes:
                    letters_with_positive_votes.append(ele)

        if len(letters_with_positive_votes) > 0:
            letters_with_negative_votes = [lett for lett in letters_present if lett not in letters_with_positive_votes
                                           and lett!='0' and lett!='-' and lett!=0]
        else:
            letters_with_negative_votes = []

        pos_counts_dict = {}
        neg_counts_dict = {}

        #count how many times positive vote letter is repeated on either side of read
        for letter in letters_with_positive_votes:
            letter_vote_multiplicity = []
            for block in alignment_blocks:
                if letter in block:
                    letter_vote_multiplicity.append(block.count(letter))
            pos_counts_dict[letter] = min(letter_vote_multiplicity)

        #count how many times negative vote letter is NOT repeated on either side of read
        for letter in letters_with_negative_votes:
            letter_vote_multiplicity = 1
            for block in alignment_blocks:
                if letter in block:
                    letter_vote_multiplicity=block.count(letter)
            neg_counts_dict[letter] = letter_vote_multiplicity

        for letter in letters_with_negative_votes:
            if letter!='*' and letter!=0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][1]+= (1 * neg_counts_dict[letter])
                else:
                    the_dict[letter][1] += (1 * neg_counts_dict[letter])

        for letter in letters_with_positive_votes:
            if letter != '*' and letter != 0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][0]+= (1* pos_counts_dict[letter])
                else:
                    the_dict[letter][0] += (1* pos_counts_dict[letter])

def update_vote_dict_only_repeats(test, the_dict):
    letters_with_positive_votes = []
    letters_with_negative_votes = []

    # for reads with no artifacts, return LRSes
    if '-' not in test:
        test_abs = [abs(ele) if isinstance(ele, int) else ele for ele in test]
        letters_present = list(set(test_abs))

        #update votes
        positive_votes =  return_l_mf_repeats(allrs(test_abs))

        if len(positive_votes) >0:

            for tup in positive_votes:
                for num in tup:
                    letters_with_positive_votes.append(num)

            letters_with_positive_votes = list(set(letters_with_positive_votes))
            letters_with_negative_votes = [lett for lett in letters_present if lett not in letters_with_positive_votes]

        else:
            letters_with_positive_votes = []
            letters_with_negative_votes = []

        #update votes

        for letter in letters_with_negative_votes:
            if letter != '*' and letter != 0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][1]+=1
                else:
                    the_dict[letter][1] += 1

        for letter in letters_with_positive_votes:
            if letter != '*' and letter != 0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][0]+=1
                else:
                    the_dict[letter][0] += 1

def update_vote_dict_only_artifacts(test, the_dict):
    letters_with_positive_votes = []
    letters_with_negative_votes = []

    #for reads with artifacts, return any abs() repeated characters across aligned sections
    if '-' in test:

        indicies_delim = [i for i, n in enumerate(test) if n=='-']
        alignment_blocks = [] #list of set of chars between asterisks
        test_abs = [abs(ele) if isinstance(ele, int) else ele for ele in test]
        letters_present = list(set(test_abs))

        if len(indicies_delim)==1:
            #handle start:
            alignment_blocks.append(list(test_abs[0: indicies_delim[0]]))

            #handle end:
            alignment_blocks.append(list(test_abs[indicies_delim[0] + 1:]))

            #do intersection of both halves////////

        elif len(indicies_delim)>1:
            # handle start sections:
            alignment_blocks.append(list(test_abs[0: indicies_delim[0]]))

            # handle middle sections
            mid_ind = indicies_delim[0:len(indicies_delim)]
            for i in range(len(mid_ind)-1):
                alignment_blocks.append(list(test_abs[mid_ind[i] + 1: mid_ind[i+1]]))

            # handle end sections
            alignment_blocks.append(list(test_abs[indicies_delim[-1] + 1:]))

            #do pairwise intersection of all sections and record, then deduplicate

        #print alignment_blocks
        intersections = []
        for i, junc_list in enumerate(alignment_blocks):
            for j, junc_list2 in enumerate(alignment_blocks):
                if i!=j:
                    intersections.append(list(set(junc_list) & set(junc_list2)))

        #print intersections
        positive_votes = [list(item) for item in set(tuple(row) for row in intersections)]
        letters_with_positive_votes = []
        for some_list in positive_votes:
            for ele in some_list:
                if ele not in letters_with_positive_votes:
                    letters_with_positive_votes.append(ele)

        #negative votes can only be cast if there is something IN positive votes
        trimmed_pos_votes = [ele for ele in letters_with_positive_votes if (ele!=0) and (ele!='*')]

        if len(trimmed_pos_votes) > 0:
            letters_with_negative_votes = [lett for lett in letters_present if lett not in trimmed_pos_votes
                                            and lett!='-' and lett!=0 and lett!='*']
        else:
            letters_with_negative_votes = []

        pos_counts_dict = {}
        neg_counts_dict = {}

        #count how many times positive vote letter is repeated on either side of read
        for letter in letters_with_positive_votes:
            letter_vote_multiplicity = []
            for block in alignment_blocks:
                if letter in block:
                    letter_vote_multiplicity.append(block.count(letter))
            pos_counts_dict[letter] = min(letter_vote_multiplicity)

        #count how many times negative vote letter is NOT repeated on either side of read
        for letter in letters_with_negative_votes:
            letter_vote_multiplicity = 1
            for block in alignment_blocks:
                if letter in block:
                    letter_vote_multiplicity=block.count(letter)
            neg_counts_dict[letter] = letter_vote_multiplicity

        for letter in letters_with_negative_votes:
            if letter!='*' and letter!=0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][1]+= (1 * neg_counts_dict[letter])
                else:
                    the_dict[letter][1] += (1 * neg_counts_dict[letter])

        for letter in letters_with_positive_votes:
            if letter != '*' and letter != 0:
                if letter not in the_dict.keys():
                    the_dict[letter] = [0,0]
                    the_dict[letter][0]+= (1* pos_counts_dict[letter])
                else:
                    the_dict[letter][0] += (1* pos_counts_dict[letter])

def strip_nontrue_junctions(alphabet_read, true_junctions):
    to_return = []
    for ele in alphabet_read:
        if ele!='*' and ele!='-' and ele!=0:
            if abs(ele) in true_junctions:
                to_return.append(ele)
        else:
            to_return.append(ele)

    return to_return

def strip_nontrue_junctions_allformats(alphabet_read, orientation_read, highlow_read, true_junctions):
    to_return_alphabet = []
    to_return_orientation = []
    to_return_highlow = []
    for alphabet, orientation, highlow in zip(alphabet_read, orientation_read, highlow_read):

        if alphabet!='*' and alphabet!='-' and alphabet!=0:
            if abs(alphabet) in true_junctions:
                to_return_alphabet.append(alphabet)
                to_return_orientation.append(orientation)
                to_return_highlow.append(highlow)
        else:
            to_return_alphabet.append(alphabet)
            to_return_orientation.append(orientation)
            to_return_highlow.append(highlow)

    return to_return_alphabet, to_return_orientation, to_return_highlow

def return_junc_transitions_absolute(stripped_read):

    junction_transitions = []

    if '-' in stripped_read:
        indicies_delim = [i for i, n in enumerate(stripped_read) if n == '-']
        alignment_blocks = []  # list of set of chars between '-' entries

        if len(indicies_delim) == 1:
            # handle start:
            alignment_blocks.append(list(stripped_read[0: indicies_delim[0]]))

            # handle end:
            alignment_blocks.append(list(stripped_read[indicies_delim[0] + 1:]))

        elif len(indicies_delim) > 1:
            # handle start sections:
            alignment_blocks.append(list(stripped_read[0: indicies_delim[0]]))

            # handle middle sections
            mid_ind = indicies_delim[0:len(indicies_delim)]
            for z in range(len(mid_ind) - 1):
                alignment_blocks.append(list(stripped_read[mid_ind[z] + 1: mid_ind[z + 1]]))

            # handle end sections
            alignment_blocks.append(list(stripped_read[indicies_delim[-1] + 1:]))

    else:
        alignment_blocks = [stripped_read]

    #print alignment_blocks

    for aln_block in alignment_blocks:
        if len(aln_block) >=2:
            for k in range(len(aln_block)-1):
                transition = aln_block[k:k+2]
                if ('*' not in transition) and (0 not in transition):
                    absolute_transition = [abs(ele) for ele in transition]
                    junction_transitions.append(absolute_transition)

    return junction_transitions

def return_junc_transitions_allformats(stripped_read_alphabet, stripped_read_orientation, stripped_read_lowhigh):

    junction_transitions_alphabet = []
    junction_transitions_orientation = []
    junction_transitions_lowhigh = []

    if '-' in stripped_read_alphabet:
        indicies_delim = [i for i, n in enumerate(stripped_read_alphabet) if n == '-']

        alignment_blocks_alphabet = []  # list of set of chars between '-' entries
        alignment_blocks_orientation = []
        alignment_blocks_lowhigh = []

        if len(indicies_delim) == 1:
            # handle start:
            alignment_blocks_alphabet.append(list(stripped_read_alphabet[0: indicies_delim[0]]))
            alignment_blocks_orientation.append(list(stripped_read_orientation[0: indicies_delim[0]]))
            alignment_blocks_lowhigh.append(list(stripped_read_lowhigh[0: indicies_delim[0]]))

            # handle end:
            alignment_blocks_alphabet.append(list(stripped_read_alphabet[indicies_delim[0] + 1:]))
            alignment_blocks_orientation.append(list(stripped_read_orientation[indicies_delim[0] + 1:]))
            alignment_blocks_lowhigh.append(list(stripped_read_lowhigh[indicies_delim[0] + 1:]))

        elif len(indicies_delim) > 1:
            # handle start sections:
            alignment_blocks_alphabet.append(list(stripped_read_alphabet[0: indicies_delim[0]]))
            alignment_blocks_orientation.append(list(stripped_read_orientation[0: indicies_delim[0]]))
            alignment_blocks_lowhigh.append(list(stripped_read_lowhigh[0: indicies_delim[0]]))

            # handle middle sections
            mid_ind = indicies_delim[0:len(indicies_delim)]
            for z in range(len(mid_ind) - 1):
                alignment_blocks_alphabet.append(list(stripped_read_alphabet[mid_ind[z] + 1: mid_ind[z + 1]]))
                alignment_blocks_orientation.append(list(stripped_read_orientation[mid_ind[z] + 1: mid_ind[z + 1]]))
                alignment_blocks_lowhigh.append(list(stripped_read_lowhigh[mid_ind[z] + 1: mid_ind[z + 1]]))

            # handle end sections
            alignment_blocks_alphabet.append(list(stripped_read_alphabet[indicies_delim[-1] + 1:]))
            alignment_blocks_orientation.append(list(stripped_read_orientation[indicies_delim[-1] + 1:]))
            alignment_blocks_lowhigh.append(list(stripped_read_lowhigh[indicies_delim[-1] + 1:]))

    else:
        alignment_blocks_alphabet = [stripped_read_alphabet]
        alignment_blocks_orientation = [stripped_read_orientation]
        alignment_blocks_lowhigh = [stripped_read_lowhigh]

    for aln_block_alphabet, aln_block_orientation, aln_block_lowhigh in \
            zip(alignment_blocks_alphabet, alignment_blocks_orientation, alignment_blocks_lowhigh):

        if len(aln_block_alphabet) >=2:
            for k in range(len(aln_block_alphabet)-1):
                transition_alphabet = aln_block_alphabet[k:k+2]
                transition_orientation = aln_block_orientation[k:k + 2]
                transition_lowhigh = aln_block_lowhigh[k:k + 2]

                if ('*' not in transition_alphabet) and (0 not in transition_alphabet):
                    #absolute_transition_alphabet = [abs(ele) for ele in transition_alphabet]
                    junction_transitions_alphabet.append(transition_alphabet)
                    junction_transitions_orientation.append(transition_orientation)
                    junction_transitions_lowhigh.append(transition_lowhigh)

    return junction_transitions_alphabet, junction_transitions_orientation, junction_transitions_lowhigh

def extract_ele_list_pairs(some_list):
    elements = []
    for pair in some_list:
        for ele in pair:
            elements.append(ele)
    return elements

def return_lcs_repeats(wholeref_read_alphabet, wholeref_read_alphabet_orientation, wholeref_read_alphabet_num_order):
    #for every read, compute l_mf_repeat between '-' entries, or whole read if no artifact present to detect repeats
    full_repeats = []
    for read, read_orientation, read_num_order in zip(wholeref_read_alphabet, wholeref_read_alphabet_orientation, wholeref_read_alphabet_num_order):

        if '-' in read:

            indicies_delim = [i for i, n in enumerate(read) if n == '-']
            alignment_blocks = []  # list of set of chars between '-' entries (LH rep)
            alignment_blocks_ori = [] # (orientation rep)
            block_locations = []
            if len(indicies_delim) == 1:
                # handle start:
                alignment_blocks.append(list(read_num_order[0: indicies_delim[0]]))
                alignment_blocks_ori.append(list(read_orientation[0: indicies_delim[0]]))
                block_locations.append([0, indicies_delim[0]])

                # handle end:
                alignment_blocks.append(list(read_num_order[indicies_delim[0] + 1:]))
                alignment_blocks_ori.append(list(read_orientation[indicies_delim[0] + 1:]))
                block_locations.append([indicies_delim[0] + 1, len(read_num_order)])

            elif len(indicies_delim) > 1:
                # handle start sections:
                alignment_blocks.append(list(read_num_order[0: indicies_delim[0]]))
                alignment_blocks_ori.append(list(read_orientation[0: indicies_delim[0]]))
                block_locations.append([0, indicies_delim[0]])

                # handle middle sections
                mid_ind = indicies_delim[0:len(indicies_delim)]
                for z in range(len(mid_ind) - 1):
                    alignment_blocks.append(list(read_num_order[mid_ind[z] + 1: mid_ind[z + 1]]))
                    alignment_blocks_ori.append(list(read_orientation[mid_ind[z] + 1: mid_ind[z + 1]]))
                    block_locations.append([mid_ind[z] + 1, mid_ind[z + 1]])

                # handle end sections
                alignment_blocks.append(list(read_num_order[indicies_delim[-1] + 1:]))
                alignment_blocks_ori.append(list(read_orientation[indicies_delim[-1] + 1:]))
                block_locations.append([indicies_delim[-1] + 1, len(read_num_order)])

        else:
            alignment_blocks = [read_num_order]
            alignment_blocks_ori = [read_orientation]
            block_locations = [[0, len(read_num_order)]]

        for j, aln_block in enumerate(alignment_blocks):
            if len(aln_block)>0:

                rs_dict, rs_loc_dict = allrs_indicies(aln_block)
                repeats, repeats_loc = return_l_mf_repeats_indicies(rs_dict, rs_loc_dict)

                rs_dict_ori, rs_loc_dict_ori = allrs_indicies(alignment_blocks_ori[j])
                repeats_ori, repeats_loc_ori = return_l_mf_repeats_indicies(rs_dict_ori, rs_loc_dict_ori)

                if len(repeats)!=0:# and repeats_loc == repeats_loc_ori:

                    #imposing limit, require two periods for uniform repeat
                    for k, repeat in enumerate(repeats):

                        if len(repeat)>=2:
                            repeat_is_legal = True
                            repeat_is_tandem = False

                            #determine if repeat contains only legal letters
                            for repeat_letter in repeat:
                                int_letter = int(repeat_letter[:-2])
                                if int_letter not in true_junctions:
                                    repeat_is_legal=False

                            #determine if this is a true "tandem" repeat
                            for l in range(len(repeats_loc[k])):
                                for m in range(len(repeats_loc[k])):
                                    if l!=m and is_overlapped(repeats_loc[k][l][0], repeats_loc[k][l][1], repeats_loc[k][m][0], repeats_loc[k][m][1]):
                                        repeat_is_tandem= True

                            if repeat_is_legal and repeat_is_tandem:

                                #print repeat, repeats_loc[k], read_num_order, read_orientation

                                orientation_block = read_orientation[block_locations[j][0]:block_locations[j][1]]
                                orientation_repeat = orientation_block[repeats_loc[k][0][0]:repeats_loc[k][0][1]]

                                num_order_block = read_num_order[block_locations[j][0]:block_locations[j][1]]
                                num_order_repeat = num_order_block[repeats_loc[k][0][0]:repeats_loc[k][0][1]]

                                alphabet_block = read[block_locations[j][0]:block_locations[j][1]]
                                alphabet_repeat = alphabet_block[repeats_loc[k][0][0]:repeats_loc[k][0][1]]

                                full_repeats.append([alphabet_repeat, orientation_repeat, num_order_repeat])

                                #print [alphabet_repeat, orientation_repeat, num_order_repeat], read
    lcs_repeat_dict = {}

    #generate repeats from full_repeats list, which contains more repeat length than necessary
    for lst in full_repeats:

        repeat_alphabet = lst[0]
        repeat_orientation = lst[1]
        repeat_lowhigh = lst[2]

        print repeat_alphabet, repeat_lowhigh, repeat_orientation, len(set(repeat_lowhigh)), len(set(repeat_orientation))

        if len(set(repeat_orientation)) <= len(set(repeat_lowhigh)):
            longest_set = len(set(repeat_lowhigh)) #max(len(set(repeat_orientation)), len(set(repeat_lowhigh)))
            repeat_alphabet_slice = repeat_alphabet[:longest_set]
            repeat_orientation_slice = repeat_orientation[:longest_set]
            repeat_lowhigh_slice = repeat_lowhigh[:longest_set]

            #print repeat_orientation
            #print repeat_orientation_slice, repeat_lowhigh_slice

            present_already = False
            #special case for single values present
            if len(repeat_lowhigh_slice) ==1:
                singlet_perm = [repeat_lowhigh_slice[0][:-2] + repeat_lowhigh_slice[0][-1] + repeat_lowhigh_slice[0][-2],
                                repeat_lowhigh_slice[0]]

                #print singlet_perm, "singlet_perm"
                for rotation in singlet_perm:
                    if tuple([rotation]) in lcs_repeat_dict.keys():
                        present_already = True
            else:
                #standard rotation
                for rotation in permute_only(repeat_lowhigh_slice):
                    if tuple(rotation) in lcs_repeat_dict.keys():
                        present_already=True
                #or standard rotation WITH LH permutation
                LH_repeat_perm = [entry[:-2] + entry[-1] + entry[-2] for entry in repeat_lowhigh_slice]
                for rotation in permute_only(LH_repeat_perm):
                    if tuple(rotation) in lcs_repeat_dict.keys():
                        present_already=True

            if not present_already:
                lcs_repeat_dict[tuple(repeat_lowhigh_slice)] = [repeat_orientation_slice, repeat_alphabet_slice]

    return lcs_repeat_dict

def return_junction_ref_boundaries(some_true_junctions):

    # sift through all junctions, storing reference alignment boundaries
    # format per key can be seen below ( [stats for low] [stats for high]
    # keys are '1' etc. for all true junctions

    junction_ref_boundaries_dict= {}

    for key in alphabet_junc_dict.keys():
        highs = []
        lows = []

        if key in some_true_junctions:
            if key not in junction_ref_boundaries_dict.keys():
                junction_ref_boundaries_dict[key] = [[],[]]

            for junction in alphabet_junc_dict[key]:
                high = max([int(junction[0]), int(junction[1])])
                low = min([int(junction[0]), int(junction[1])])

                highs.append(high)
                lows.append(low)

            #print highs
            km_high = KMeans(n_clusters=2)
            km_low = KMeans(n_clusters=2)
            km_high.fit(np.array(highs).reshape(-1,1))
            km_low.fit(np.array(lows).reshape(-1,1))

            high_cluster_centers = km_high.cluster_centers_
            low_cluster_centers = km_low.cluster_centers_

            high_labels = km_high.labels_
            low_labels = km_low.labels_

            high_clusters = []
            low_clusters = []

            high_cluster_labels = []
            low_cluster_labels = []

            for i in set(high_labels):
                high_class_mask = (high_labels ==i)
                high_clusters.append(list(np.asarray(highs)[high_class_mask]))
                high_cluster_labels.append(i)

            for j in set(low_labels):
                low_class_mask = (high_labels ==i)
                low_clusters.append(list(np.asarray(lows)[low_class_mask]))
                low_cluster_labels.append(j)

            high_high_label = high_cluster_labels[np.argmax(np.array([np.mean(lst) for lst in high_clusters]))]
            high_low_label = high_cluster_labels[np.argmin(np.array([np.mean(lst) for lst in high_clusters]))]

            low_high_label = low_cluster_labels[np.argmax(np.array([np.mean(lst) for lst in low_clusters]))]
            low_low_label = low_cluster_labels[np.argmin(np.array([np.mean(lst) for lst in low_clusters]))]

            high_high_mean = max([ele[0] for ele in high_cluster_centers])
            high_low_mean = min([ele[0] for ele in high_cluster_centers])

            low_high_mean = max([ele[0] for ele in low_cluster_centers])
            low_low_mean = min([ele[0] for ele in low_cluster_centers])

            high_low_std = np.std(np.array(highs)[(high_labels == high_low_label)])
            high_high_std = np.std(np.array(highs)[(high_labels == high_high_label)])

            low_low_std = np.std(np.array(lows)[(low_labels == low_low_label)])
            low_high_std = np.std(np.array(lows)[(low_labels == low_high_label)])

            #print high_high_mean, high_high_std

            #, np.std(np.asarray(lows))
            # , np.std(np.asarray(highs))
            junction_ref_boundaries_dict[key][0].extend([int(low_low_mean), int(low_high_mean), int(np.mean(np.asarray(lows)))])
            junction_ref_boundaries_dict[key][1].extend([int(high_low_mean), int(high_high_mean), int(np.mean(np.asarray(highs)))])

            #now append standard deviations, so we have low stats, high stats, low stat stdev, high stat stdev
            junction_ref_boundaries_dict[key].append([low_low_std, low_high_std, np.std(np.asarray(lows))])
            junction_ref_boundaries_dict[key].append([high_low_std, high_high_std, np.std(np.asarray(highs))])

    return junction_ref_boundaries_dict

def return_junction_inversion_status(some_true_junctions):

    # sift through all examples of each junction flavour and record whether or not it is inverted
    junctions_inverted= {}

    for key in alphabet_junc_dict.keys():

        if key in some_true_junctions:

            junctions_inverted[key] = alphabet_junc_dict[key][0][2]

    return junctions_inverted

def add_return_ref_locations_dict(junction_ref_boundaries_dict, input_dict):

    for key in input_dict.keys():

        #print key
        #create a template key to make constructing repeat easier (adding first entry to end of existing key)
        template_key = [ele for ele in key]
        template_key.append(key[0])

        template_key_sign = [ele for ele in input_dict[key][0]]
        template_key_sign.append(input_dict[key][0][0])

        #take all interior transitions in templates to construct alignments using junction_ref_boundaries_dict
        input_dict[key].append([])

        for i in range(len(template_key)-1):
            curr_entry = template_key[i]
            next_entry = template_key[i+1]

            curr_entry_num = int(curr_entry[:-2])
            next_entry_num = int(next_entry[:-2])

            curr_entry_level = curr_entry[-1]
            next_entry_level = next_entry[-2]

            #define orientation of alignment
            if template_key_sign[i][-1] != template_key_sign[i+1][-2]:
                print "alignment orientation error", template_key
            else:
                alignment_sign = template_key_sign[i][-1]

            #consider what entries to use based on level transition
            if curr_entry_level == 'L':
                curr_ref_index = 0
            else:
                curr_ref_index = 1

            if next_entry_level == 'L':
                next_ref_index = 0
            else:
                next_ref_index = 1

            if (curr_entry_level == next_entry_level) and (curr_entry_num == next_entry_num) and alignment_sign=='-':
                start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][1]
                end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][0]

                start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index+2][1]
                end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index+2][0]

            elif (curr_entry_level == next_entry_level) and (curr_entry_num == next_entry_num) and alignment_sign=='+':
                start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][0]
                end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][1]

                start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index+2][0]
                end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index+2][1]

            else:
                start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][2]
                end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][2]

                start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index+2][2]
                end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index+2][2]

            # if [start,end, start_std, end_std] not in input_dict[key][-1]:
            input_dict[key][-1].append([start,end, start_std, end_std])

            #     #being here would indicate a noncontiguous repeat, so remove it...
            # else:
            #     input_dict.pop(key)

    return input_dict

#for single reads, not detected repeats
def add_return_ref_locations_dict_example(junction_ref_boundaries_dict, input_dict):
    for key in input_dict.keys():

        #create a template key to make constructing repeat easier (adding first entry to end of existing key)
        template_key = [ele for ele in key]

        template_key_sign = [ele for ele in input_dict[key][0]]

        #take all interior transitions in templates to construct alignments using junction_ref_boundaries_dict
        input_dict[key].append([])

        for i in range(len(template_key)-1):
            curr_entry = template_key[i]
            next_entry = template_key[i+1]

            curr_entry_num = int(curr_entry[:-2])
            next_entry_num = int(next_entry[:-2])

            curr_entry_level = curr_entry[-1]
            next_entry_level = next_entry[-2]

            #define orientation of alignment
            if template_key_sign[i][-1] != template_key_sign[i+1][-2]:
                print "alignment orientation error"
            else:
                alignment_sign = template_key_sign[i][-1]

            #consider what entries to use based on level transition
            if curr_entry_level == 'L':
                curr_ref_index = 0
            else:
                curr_ref_index = 1

            if next_entry_level == 'L':
                next_ref_index = 0
            else:
                next_ref_index = 1

            if (curr_entry_level == next_entry_level) and (curr_entry_num == next_entry_num) and alignment_sign=='-':
                start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][1]
                end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][0]

                start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index+2][1]
                end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index+2][0]

            elif (curr_entry_level == next_entry_level) and (curr_entry_num == next_entry_num) and alignment_sign=='+':
                start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][0]
                end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][1]

                start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index+2][0]
                end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index+2][1]

            else:
                start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][2]
                end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][2]

                start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index+2][2]
                end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index+2][2]


            input_dict[key][-1].append([start,end, start_std, end_std])

    return input_dict

#iterate through all reads, and grab those that correspond to primary junction transitions to extract mixed alignments
def return_mixed_alignment_count_dict(wholeref_read_alphabet, wholeref_read_alphabet_orientation, wholeref_read_alphabet_num_order, primary_junctions_trans):
    if primary_mixed_repeat==True:
        mixed_alignment_count_dict = {}
        for read_alphabet, read_orientation, read_low_high in zip(wholeref_read_alphabet, wholeref_read_alphabet_orientation, wholeref_read_alphabet_num_order):
            stripped_alphabet, stripped_orientation, stripped_lowhigh = strip_nontrue_junctions_allformats(read_alphabet, read_orientation, read_low_high, true_junctions)
            transitions_alphabet, transitions_orientation, transitions_lowhigh = return_junc_transitions_allformats(stripped_alphabet, stripped_orientation, stripped_lowhigh)

            if len(transitions_alphabet)>0:
                for i, trans in enumerate(transitions_alphabet):
                    if [abs(junc) for junc in trans] in primary_junctions_trans:

                        aln_template = transitions_lowhigh[i]
                        aln_template_sign = transitions_orientation[i]

                        curr_entry = aln_template[0]
                        next_entry = aln_template[1]

                        curr_entry_num = int(curr_entry[:-2])
                        next_entry_num = int(next_entry[:-2])

                        curr_entry_level = curr_entry[-1]
                        next_entry_level = next_entry[-2]

                        # define orientation of alignment
                        if aln_template_sign[0][-1] != aln_template_sign[1][-2]:
                            print
                            "alignment orientation error"
                        else:
                            alignment_sign = aln_template_sign[0][-1]

                        # consider what entries to use based on level transition
                        if curr_entry_level == 'L':
                            curr_ref_index = 0
                        else:
                            curr_ref_index = 1

                        if next_entry_level == 'L':
                            next_ref_index = 0
                        else:
                            next_ref_index = 1

                        if curr_entry_level == next_entry_level and (curr_entry_num == next_entry_num) and alignment_sign == '-':
                            start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][1]
                            end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][0]

                            start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index + 2][1]
                            end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index + 2][0]

                        elif curr_entry_level == next_entry_level and (curr_entry_num == next_entry_num) and  alignment_sign == '+':
                            start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][0]
                            end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][1]

                            start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index + 2][0]
                            end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index + 2][1]

                        else:
                            start = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index][2]
                            end = junction_ref_boundaries_dict[next_entry_num][next_ref_index][2]

                            start_std = junction_ref_boundaries_dict[curr_entry_num][curr_ref_index + 2][2]
                            end_std = junction_ref_boundaries_dict[next_entry_num][next_ref_index + 2][2]

                            #print aln_template, aln_template_sign, start, end, start_std, end_std

                        potential_key = [aln_template[0][:-2] + aln_template[0][-1], aln_template[1][:-2] + aln_template[1][-2]]

                        in_dict=False
                        to_append_key = ""
                        for key in mixed_alignment_count_dict.keys():

                            if potential_key in permute_only(list(key)):
                                in_dict=True
                                to_append_key=key

                        if not in_dict:
                            mixed_alignment_count_dict[tuple(potential_key)] = [[start, end, start_std, end_std],1]

                        else:
                            mixed_alignment_count_dict[to_append_key][-1]+=1
        return mixed_alignment_count_dict

def juncs_to_vect(junctions, true_junction_in_repeats):
    length = len(true_junction_in_repeats)
    vect = np.zeros(length)
    for junc in junctions:
        index = true_junction_in_repeats.index(junc)
        vect[index]=1

    return list(vect)

# given a read in alphabet format and its endpoints, extract alignments in list format
def get_alignments_from_read(read_as_alphabet, endpoints, true_junctions, read_juncformat):
    alignments_to_return = []
    indicies_to_remove = []
    read_as_alphabet_abs = []

    #populate read_as_alphabet_abs
    for junc in read_as_alphabet:
        if isinstance(junc, int):
            read_as_alphabet_abs.append(abs(junc))
        else:
            read_as_alphabet_abs.append(junc)

    #record locations where junctions are not from true_junctions
    for i, junc in enumerate(read_as_alphabet_abs):
        if isinstance(junc, int):
            if junc not in true_junctions:
                indicies_to_remove.append(i)
        else:
            if junc!='-':
                indicies_to_remove.append(i)

    #remove junctions in juncformat that we don't believe are real (spurious insertion)
    cleaned_juncformat = []
    cleaned_alphabetformat = []
    for j, junc in enumerate(read_juncformat):
        if j not in indicies_to_remove:
            cleaned_juncformat.append(read_juncformat[j])
            cleaned_alphabetformat.append(read_as_alphabet_abs[j])

    #handle interior alignments
    if any(ele in cleaned_alphabetformat for ele in true_junctions):
        if cleaned_juncformat>=2:
            for i in range(len(cleaned_juncformat)-1):
                alignments_to_return.append([cleaned_juncformat[i][6], cleaned_juncformat[i+1][5]])

        #handle start
        if 0 not in indicies_to_remove:
            alignments_to_return.append([endpoints[0],cleaned_juncformat[0][5]])

        #handle end
        if len(read_as_alphabet)-1 not in indicies_to_remove:
            alignments_to_return.append([cleaned_juncformat[-1][6], endpoints[1]])

        return alignments_to_return

def return_fractional_nonoverlap(query_list, reference_list):
    nonoverlap_length = []
    reference_span = 0

    #perform pairwise query to reference alignment test
    for q_list in query_list:
        per_q_nonoverlap = []
        for r_list in reference_list:
            q_range = set(range(min(q_list[0], q_list[1]), max(q_list[0], q_list[1]) +1))
            r_range = set(range(min(r_list[0], r_list[1]), max(r_list[0], r_list[1]) + 1))

            symm_diff = (q_range ^ r_range) & q_range
            per_q_nonoverlap.append(len(symm_diff))

        nonoverlap_length.append(min(per_q_nonoverlap))

    #now sum span of reference ranges
    for r_list in reference_list:
        reference_span+= abs(max(r_list[0], r_list[1]) - min(r_list[0], r_list[1]))

    nonoverlap_frac = sum(nonoverlap_length)/float(reference_span)

    return nonoverlap_frac

complimentary={'A':'T', 'C': 'G', 'T':'A', 'G':'C'}

'''Reading in bam files...'''
bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
#bamfile2 = pysam.AlignmentFile(sys.argv[2], "rb")
print sys.argv[1].strip().split('/')[-1][0:-4]
output_name = sys.argv[1].strip().split('/')[-1][0:-4].split('_')[0] + "_" + sys.argv[1].strip().split('/')[-1][0:-4].split('_')[1]
print output_name
complete_read_list = []  # all entries in BAM, included supplementary and secondary reads
unique_read_list = []  # only primary reads in BAM

'''Reading in invdup read locations aggregrated across all ChrI-ChrIII reads'''
null_invdup_locations = []
binwidth=0.025
with open('inverted_duplication_read_locations_chrI-III.txt', 'r') as file_in:
    for line in file_in:
        null_invdup_locations.append(float(line.strip().split()[0]))

prob, invdup_readloc_bins = np.histogram(null_invdup_locations, bins=np.arange(0, 1+ binwidth, binwidth), density=True)
invdup_loc_prob = binwidth*prob #if invdup_loc_prob <5% for a particular region, include on potential list of real invdups

# construct read lists
for entry in bamfile:
    lor = entry.query_length
    if lor > 0:
        if entry.is_unmapped == False and entry.is_supplementary == False and entry.is_secondary == False:
            #secondary set of criteria, if primary is not from a clipped region
            if not "offset" in entry.query_name.strip().split('_'):
                unique_read_list.append([entry, entry.query_length])
        if entry.is_unmapped == False and entry.is_secondary == False:
            complete_read_list.append([entry, entry.query_length])

total = len(unique_read_list)  # number of unique reads
tail_read_count = total - int(0.0 * total)  # how many unique read names we want in the tail of the size distribution

# list of unique primary reads
sorted_unique_read_list = sorted(unique_read_list, key=lambda x: x[1], reverse=True)
tail_primary_reads = sorted_unique_read_list[:tail_read_count]

returned_junctions_wholeref = []
returned_read_endpoints_wholeref = []
returned_read_juncformat_wholeref = []
returned_read_lengths_wholeref = []
invdup_locations_wholeref = []
invdup_fractional_locations_wholeref = [] #all reads (including those with multiple invdup signals)
invdups_reflocations_wholeref = []
invdup_locations_multiples_wholeref = [] #derived from reads with multiple invdup signals
invdup_fractional_locations_multiples_wholeref = []
invdups_reflocations_multiples_wholeref = []

for read_p in tail_primary_reads:
    junction_in_primary = []
    subalignments = []
    read_length = read_p[0].infer_read_length()
    primary_name = read_p[0].query_name
    #print primary_name

    #extract primary sequence
    if read_p[0].is_reverse:
        sequence = read_p[0].query_sequence
        sequence = complimentary_strand(sequence)
    else:
        sequence = read_p[0].query_sequence

    for read in complete_read_list:
        #if a non clipped segment
        if read_p[0].query_name == read[0].query_name and read[0].mapping_quality>=20:
            read_name = read_p[0].query_name
            clips = []
            clip_positions = []
            cigar_tuples = read[0].cigartuples
            blocks = read[0].get_blocks()
            matches_insertions_wholeread = []

            if read[0].is_reverse:
                cigar_tuples = cigar_tuples[::-1]
                blocks = blocks[::-1]

            for i, tup in enumerate(cigar_tuples):
                # if either soft/hard clipped (handles primary read which is only soft-clipped)
                if tup[0] == 4 or tup[0] == 5:
                    clips.append(int(tup[1]))
                    clip_positions.append(i)

            if len(clips)==2:
                qstart = clips[0] #starting position of segment in read coordinates
                qend = read_length - clips[-1]
            elif len(clips)==1:
                if clip_positions[0]==0:
                    qstart = clips[0]  # starting position of segment in read coordinates
                    qend = read_length
                else:
                    qstart = 0 # starting position of segment in read coordinates
                    qend = read_length - clips[-1]
            else:
                qstart = 0  # starting position of segment in read coordinates
                qend = read_length

            ref_start = read[0].reference_start
            ref_end = read[0].reference_end

            if read[0].is_reverse:
                subalignments.append([qstart, qend, ref_end, ref_start, read[0].is_reverse, read[0].mapq, read_name])
            else:
                subalignments.append([qstart, qend, ref_start, ref_end, read[0].is_reverse, read[0].mapq, read_name])

        #if this is from a clipped region, make sure to include offsets
        elif read_p[0].query_name == read[0].query_name.split("offset")[0][:-1] and read[0].mapping_quality>=20:
            read_name = read_p[0].query_name
            clips = []
            clip_positions = []
            cigar_tuples = read[0].cigartuples
            blocks = read[0].get_blocks()
            matches_insertions_wholeread = []

            start_offset = int(read[0].query_name.split("_")[-2])
            end_offset = int(read[0].query_name.split("_")[-1])

            if read[0].is_reverse:
                cigar_tuples = cigar_tuples[::-1]
                blocks = blocks[::-1]

            for i, tup in enumerate(cigar_tuples):
                # if either soft/hard clipped (handles primary read which is only soft-clipped)
                if tup[0] == 4 or tup[0] == 5:
                    clips.append(int(tup[1]))
                    clip_positions.append(i)

            if len(clips) == 2:

                qstart = clips[0]  + start_offset# starting position of segment in read coordinates
                qend = end_offset - clips[-1]

            elif len(clips) == 1:

                if clip_positions[0] == 0:
                    qstart = clips[0]  + start_offset # starting position of segment in read coordinates
                    qend = end_offset
                else:
                    qstart = start_offset # starting position of segment in read coordinates
                    qend = end_offset - clips[-1]
            else:
                qstart = start_offset # starting position of segment in read coordinates
                qend = end_offset

            ref_start = read[0].reference_start
            ref_end = read[0].reference_end

            if read[0].is_reverse:
                subalignments.append([qstart, qend, ref_end, ref_start, read[0].is_reverse, read[0].mapq, read_name])
            else:
                subalignments.append([qstart, qend, ref_start, ref_end, read[0].is_reverse, read[0].mapq, read_name])

    #sort by starting position in read
    sorted_subalignments = sorted(subalignments, key=lambda x:x[0])

    #cull subalignments if less than 300bp (likely spurious) because they don't contain enough content
    to_remove = []
    for alignment in sorted_subalignments:
        if abs(alignment[0]-alignment[1])< min_alignment_length:
            to_remove.append(alignment)

    sorted_subalignments = [ele for ele in sorted_subalignments if ele not in to_remove]

    #extending list that contains junctions across all reads for one strain

    if sorted_subalignments != []:

        junctions_reported, junction_locations_reported, endpoints = return_junctions_and_locations(sorted_subalignments, del_ins_threshold)

        returned_junctions_wholeref.extend(junctions_reported)

        returned_read_juncformat_wholeref.append(junctions_reported)
        returned_read_endpoints_wholeref.append(endpoints)

        returned_read_lengths_wholeref.append(read_length)


        #now begin compiling lists of invdups, for all reads and separately for reads that contain multiple invdups
        invdup_locs_per_read = []
        invdup_reflocs_per_read = []

        if len(junctions_reported) > 0:
            contains_inversion = False
            for junc, junc_loc in zip(junctions_reported, junction_locations_reported):
                if junc[2]==True and abs(junc[0] - junc[1]) < 1000: #our criteria for invdup
                    invdup_locs_per_read.append(junc_loc)
                    invdup_reflocs_per_read.append(junc)

                    invdup_locations_wholeref.append(junc_loc) #store true location in read
                    invdup_fractional_locations_wholeref.append(junc_loc/float(read_length))#store fractional location in read
                    invdups_reflocations_wholeref.append(junc)

                    contains_inversion = True

            # if we have multiple invdups per read
            if contains_inversion == True and len(invdup_locs_per_read) >=1:
                #print invdup_reflocs_per_read, primary_name, read_length
                for junc_loc in invdup_locs_per_read:
                    invdup_locations_multiples_wholeref.append(junc_loc)
                    invdup_fractional_locations_multiples_wholeref.append(junc_loc/float(read_length))
                for junc in invdup_reflocs_per_read:
                    invdups_reflocations_multiples_wholeref.append(junc)

#process wholeref lists, removing invdups based on fractional location that have prob(artifact) >0.10
invdups_post_locfilter = []
artifacts_invdups_post_locfilter = []
artifacts_invdups_post_locfilter_intersection = []
for junc, junc_loc in zip(invdups_reflocations_wholeref, invdup_fractional_locations_wholeref):
    if return_prob_invdup_location_artifact(junc_loc)<0.01: #so this number is precisely the false positive rate for junctions to be called REAL invdups if separation less than 1kbp
        invdups_post_locfilter.append(junc)
    else:
        artifacts_invdups_post_locfilter.append(junc)

#print invdups_reflocations_multiples_wholeref

#calculate intersection of list of all invdups that pass filtering and reads containing multiple invdups
invdups_post_locfilter_intersection = [list(x) for x in set(tuple(x) for x in invdups_reflocations_multiples_wholeref).intersection(set(tuple(x) for x in invdups_post_locfilter))]

artifacts_invdups_post_intersection = [junc for junc in invdups_post_locfilter if junc not in invdups_post_locfilter_intersection]

artifacts_invdups_post_locfilter.extend(artifacts_invdups_post_intersection)
artifacts_invdups_post_locfilter_intersection = artifacts_invdups_post_locfilter

currDir = os.path.dirname(os.path.realpath(__file__))
path = currDir +os.sep

#now apply clustering just on invdups, which acts as a third round of invdup filtering (first two select against artifacts, so this relies on a lower frequency of artifacts as background noise)
invdups_post_locfilter_intersection_clustering = []
artifacts_invdups_post_locfilter_intersection_clustering = []
if len(invdups_post_locfilter_intersection) >= 1:
    invdup_eps = min(compute_optimal_eps(invdups_post_locfilter_intersection), 1000)#1000
    print "eps", invdup_eps
    invdup_labels, invdup_core_samples_mask = compute_DBSCAN(invdup_eps, 3, invdups_post_locfilter_intersection)

    for k in set(invdup_labels):
        class_member_mask = (invdup_labels == k)
        if k!=-1:
            #print class_member_mask
            #print invdups_post_locfilter_intersection
            invdups_post_locfilter_intersection_clustering.extend(np.array(invdups_post_locfilter_intersection)[class_member_mask].tolist())
        else:
            artifacts_invdups_post_locfilter_intersection.extend(np.array(invdups_post_locfilter_intersection)[class_member_mask].tolist())

artifacts_invdups_post_locfilter_intersection_clustering = artifacts_invdups_post_locfilter_intersection

#print invdups_post_locfilter_intersection_clustering
#print len(invdups_post_locfilter_intersection_clustering)

#add completely filtered invdups output to list containing inversions beyond our invdup definition of <1000 separation
Z = [junc for junc in returned_junctions_wholeref if (junc[2]==True and abs(junc[1] - junc[0]) >1000)]
Z.extend(invdups_post_locfilter_intersection_clustering)
Y = [junc for junc in returned_junctions_wholeref if (junc[2]==False)] #normal junction

'''Start by clustering each class type'''
all_class_counts = []
#for inversions class
if len(Z) >= 1:
    inv_eps = cluster_merging_eps#1000 #min(compute_optimal_eps(Z), 1000)#1000
    print "eps", inv_eps
    inv_labels, inv_core_samples_mask = compute_DBSCAN(inv_eps, 4, Z)

#for normal junctions
if len(Y) >= 1:
    print compute_optimal_eps(Y)
    normal_eps = cluster_merging_eps #1000 #min(1000, compute_optimal_eps(Y))
    normal_labels, normal_core_samples_mask = compute_DBSCAN(normal_eps, 4, Y)

'''Next filter out clusters that don't reach read count contribution threshold'''

#for inverted junctions
inv_labels_readcount, inv_core_samples_mask_readcount = filter_cluster_read_support(inv_labels, Z, inv_core_samples_mask, min_read_threshold)

#for normal junctions
normal_labels_readcount, normal_core_samples_mask_readcount = filter_cluster_read_support(normal_labels, Y, normal_core_samples_mask, min_read_threshold)

'''Do cluster size counting for each type'''

#for inverted junctions
all_class_counts.extend(count_cluster_sizes(inv_labels_readcount, inv_core_samples_mask_readcount, Z))
n_clusters_inverted = len(set(inv_labels_readcount)) - (1 if -1 in inv_labels_readcount else 0)

#for normal junctions
all_class_counts.extend(count_cluster_sizes(normal_labels_readcount, normal_core_samples_mask_readcount, Y))
n_clusters_normal = len(set(normal_labels_readcount)) - (1 if -1 in normal_labels_readcount else 0)

'''Assign unique integer to clusters that pass all filtering thresholds'''
alphabet_junc_dict = {}
running_cluster_count = 0
#first populate with normal junctions
unique_normal_junctions = set(normal_labels_readcount)

#print all_class_counts

for k in unique_normal_junctions:
    junction_array = np.array(Y)
    if k!=-1:
        class_member_mask = (normal_labels_readcount == k)
        cluster = junction_array[class_member_mask]

        if running_cluster_count==0:
            running_cluster_count = 1
            alphabet_junc_dict[running_cluster_count] = cluster
            running_cluster_count+=1

        else:
            alphabet_junc_dict[running_cluster_count] = cluster
            running_cluster_count+=1

#next populate with inverted junctions
unique_inverted_junctions = set(inv_labels_readcount)

for k in unique_inverted_junctions:
    junction_array = np.array(Z)
    if k!=-1:
        class_member_mask = (inv_labels_readcount == k)
        cluster = junction_array[class_member_mask]

        alphabet_junc_dict[running_cluster_count] = cluster
        running_cluster_count+=1


wholeref_read_alphabet = []
wholeref_read_name_alphabet = []
wholeref_read_alphabet_orientation = []
wholeref_read_alphabet_num_order = []
wholeref_read_endpoints = []

dict_cluster_voting = {} #stores [#num positive votes, #num negative votes] for each cluster
dict_cluster_voting_artifacts = {}
dict_cluster_voting_repeats = {}
dict_orientation_string = {False:'+',True:'-'} #stores encoding for strand

printed=False
for read, endpoints in zip(returned_read_juncformat_wholeref, returned_read_endpoints_wholeref):
    read_in_alphabet = []
    read_in_alphabet_orientation = [] #orientation of junction, +1 False (pos strand), -1 True (neg strand)
    read_in_alphabet_num_order = [] #order of numbers (+1 high, -1 low) if enough discrepancy between the two

    for junc in read:

        orientation_transition = [junc[3], junc[4]]
        refloc_transition = [junc[5], junc[6]]
        identity = next((k for k,v in alphabet_junc_dict.items() if (np.asarray(junc, dtype='S36')==v).all(1).any()),0)
        identity_alternate = next((k for k,v in alphabet_junc_dict.items() if (np.asarray(junc, dtype='S36')==v).all(1).any()),0)
        #for inverted juncs
        if orientation_transition[0]!=orientation_transition[1]:

            #(false, low), (true, high) or (false, high), (true, low)
            #(false, false), (true, true), or  (false, true), (true, false)
            #[False, True] and [call Low False] we call positive for inverted junctions
            if refloc_transition[0] < refloc_transition[1]:
                association = [orientation_transition[0], False]
            elif refloc_transition[0] > refloc_transition[1]:
                association = [orientation_transition[0], True]
            else:
                association = [False, False]

            if association[0] == association[1]:
                identity = identity
            else:
                identity = -1* identity

        #for noninverted junctions, set false false to +1, true true to -1
        elif orientation_transition[0] == False:
            identity = identity

        else:
            identity = -1 * identity  ####what is going on here?

        #specifically identifying artifact invdup for all alphabet schemes
        if junc in artifacts_invdups_post_locfilter_intersection_clustering:
            identity = '-'
            identity_orientation = '-'
            identity_num_order = '-'
        elif identity==0 and junc[3]!=junc[4]:
            identity = '*'
            identity_orientation = '*'
            identity_num_order = '*'

        '''Record alternate alphabets that will be used in repeat reconstruction'''

        identity_orientation = str(identity_alternate) + dict_orientation_string[orientation_transition[0]] + dict_orientation_string[orientation_transition[1]]

        if refloc_transition[0] < refloc_transition[1]:
            identity_num_order = str(identity_alternate) + 'L' + 'H'

        else:
            identity_num_order = str(identity_alternate) + 'H' + 'L'

        #if junction is already in read, compare orientation to what exists
        read_in_alphabet.append(identity)
        read_in_alphabet_orientation.append(identity_orientation)
        read_in_alphabet_num_order.append(identity_num_order)

    wholeref_read_alphabet.append(read_in_alphabet)
    wholeref_read_alphabet_orientation.append(read_in_alphabet_orientation)
    wholeref_read_alphabet_num_order.append(read_in_alphabet_num_order)
    wholeref_read_endpoints.append([endpoints[0][0], endpoints[1][0]])

    wholeref_read_name_alphabet.append(junc[-1])

    # if len(read_in_alphabet)>=4:
    #     #print read_in_alphabet, junc[-1], read
    #     print read_in_alphabet, read_in_alphabet_orientation, read_in_alphabet_num_order, endpoints, read

    '''DO VOTING SCHEME, which updates voting dict'''
    update_vote_dict(read_in_alphabet, dict_cluster_voting)
    update_vote_dict_only_artifacts(read_in_alphabet, dict_cluster_voting_artifacts)
    update_vote_dict_only_repeats(read_in_alphabet, dict_cluster_voting_repeats)

# print dict_cluster_voting_artifacts
# print dict_cluster_voting_repeats

#for all dict keys, record ones where it passes artifact test
true_junctions = []
for key in dict_cluster_voting_artifacts.keys():
    if dict_cluster_voting_artifacts[key][0] > dict_cluster_voting_artifacts[key][1]:
        true_junctions.append(key)

#if key isn't in artifacts test results at all, then ask whether it passes repeat test
for key in dict_cluster_voting_repeats.keys():
    if key not in dict_cluster_voting_artifacts.keys():
        if dict_cluster_voting_repeats[key][0] > dict_cluster_voting_repeats[key][1]:
            true_junctions.append(key)

junction_counts_dict = {}
for key in alphabet_junc_dict.keys():
    if key in true_junctions:
        junction_counts_dict[key]=len(alphabet_junc_dict[key])

axis_length = len(true_junctions)
transition_matrix = np.zeros((axis_length, axis_length))

true_junctions_sorted = sorted(true_junctions)
junc_to_index_dict = {}

#print true_junctions_sorted

#construct index used to parse junc transition matrix
for i in range(len(true_junctions_sorted)):
    #first do negative versions of junctions
    junc_to_index_dict[true_junctions_sorted[i]] = i

#now populate the transition matrix using this index
for read in wholeref_read_alphabet:
    junctions = return_junc_transitions_absolute(strip_nontrue_junctions(read, true_junctions))
    for junction in junctions:
        if len(junction)==2:
            transition_matrix[junc_to_index_dict[junction[0]]][junc_to_index_dict[junction[1]]]+=1

# print transition_matrix
#print transition_matrix[junc_to_index_dict[11]][junc_to_index_dict[12]]

#parse matrix and do transition merging based on perm and binom criteria
trans_merged_binom_share= [] #master list containing merged transitions based permutation and binom test
count_merged_binom_share = []
for r in range(transition_matrix.shape[0]):
    for c in range(transition_matrix.shape[1]):
        curr_entry_count = transition_matrix[r][c]
        curr_entry_trans = [[key for key, value in junc_to_index_dict.items() if value==r][0],
                            [key for key, value in junc_to_index_dict.items() if value==c][0]]

        for r2 in range(transition_matrix.shape[0]):
            for c2 in range(transition_matrix.shape[1]):

                if r!=r2 and c!=c2:
                    next_entry_count = transition_matrix[r2][c2]
                    next_entry_trans = [[key for key, value in junc_to_index_dict.items() if value==r2][0],
                                        [key for key, value in junc_to_index_dict.items() if value==c2][0]]

                    pairwise_count = next_entry_count + curr_entry_count
                    if pairwise_count !=0:
                        bin_frac_off_u = abs(0.5* pairwise_count - curr_entry_count)/float(0.5*pairwise_count)
                        share_junc = bool(set(curr_entry_trans) & set(next_entry_trans))

                        if share_junc and bin_frac_off_u < binomial_merging_frac:

                            if len(trans_merged_binom_share) ==0:
                                trans_merged_binom_share.append([next_entry_trans, curr_entry_trans])
                                count_merged_binom_share.append([next_entry_count, curr_entry_count])

                            else:
                                add_new=True

                                #if pairs share a junction, AND pass binom test then merge into existing lists
                                for count, lst in zip(count_merged_binom_share, trans_merged_binom_share):

                                    for i in range(len(count)):
                                        share_junc_curr = bool(set(curr_entry_trans) & set(lst[i]))
                                        share_junc_next = bool(set(next_entry_trans) & set(lst[i]))

                                        bin_frac_next = abs(0.5 * (count[i] + next_entry_count) - count[i]) / float(
                                            0.5 * (count[i] + next_entry_count))
                                        bin_frac_curr = abs(0.5 * (count[i] + curr_entry_count) - count[i])/float(
                                            0.5 *(count[i] + curr_entry_count))

                                        #print share_junc_curr, share_junc_next, bin_frac_next, bin_frac_curr

                                        if share_junc_curr and bin_frac_curr < binomial_merging_frac:
                                            if next_entry_trans in lst and curr_entry_trans not in lst:
                                                lst.extend([curr_entry_trans])
                                                count.extend([curr_entry_count])

                                            elif curr_entry_trans in lst and next_entry_trans not in lst:
                                                lst.extend([next_entry_trans])
                                                count.extend([next_entry_count])

                                            add_new = False

                                        elif share_junc_next and bin_frac_next < binomial_merging_frac:
                                            if next_entry_trans in lst and curr_entry_trans not in lst:
                                                lst.extend([curr_entry_trans])
                                                count.extend([curr_entry_count])

                                            elif curr_entry_trans in lst and next_entry_trans not in lst:
                                                lst.extend([next_entry_trans])
                                                count.extend([next_entry_count])

                                            add_new = False

                                #if pairs pass binom and share entry, but don't merge with existing cluster, add new
                                if add_new:
                                    trans_merged_binom_share.append([next_entry_trans, curr_entry_trans])
                                    count_merged_binom_share.append([next_entry_count, curr_entry_count])


#output of trans_merged_binom_share can be repeated if repeats contain more than 2 junctions, so
#remove duplicates even if order is different
lst_to_remove = []
counts_to_remove = []
cleaned_trans_associations = []
cleaned_trans_counts = []
for lst, counts in zip(trans_merged_binom_share, count_merged_binom_share):
    for lst2 in trans_merged_binom_share:
        if lst!=lst2:
            shared_count = 0
            for ele in lst2:
                if ele in lst:
                    shared_count+=1
            if shared_count == len(lst) and lst2 not in lst_to_remove:
                lst_to_remove.append(lst2)
                counts_to_remove.append(counts)

for lst, counts in zip(trans_merged_binom_share, count_merged_binom_share):
    if lst not in lst_to_remove:
        cleaned_trans_associations.append(lst)
        cleaned_trans_counts.append(counts)

#now take maximum in junction encoding to represent "primary junctions", unless any number of
# entries in matrix have higher than average count
if len(cleaned_trans_associations)!=0:
    cleaned_trans_associations, cleaned_trans_counts = (list(t) for t in  zip(*sorted(
        zip(cleaned_trans_associations, cleaned_trans_counts), key=lambda pair: sum(pair[1]), reverse=True )))

    primary_junctions_trans = cleaned_trans_associations[0]

    #check if the out output constitutes a repeat, if not, discard it so that singleton is only option
    if len(set(extract_ele_list_pairs(primary_junctions_trans))) == 2 and len(primary_junctions_trans) == 2:
        primary_junctions_trans = primary_junctions_trans

    elif len(set(extract_ele_list_pairs(primary_junctions_trans))) == 1 and len(primary_junctions_trans) == 1:
        primary_junctions_trans = primary_junctions_trans

    elif 2* len(set(extract_ele_list_pairs(primary_junctions_trans))) == len(primary_junctions_trans):
        primary_junctions_trans = primary_junctions_trans

    else:
        primary_junctions_trans = []

else:
    primary_junctions_trans = []
    cleaned_trans_counts = [[0],[0]]
    cleaned_trans_associations = [[],[]]

#print cleaned_trans_associations, "clean here first"

# note that we can infer other potential repeats here...which should be used to corroborate complex LCS repeats
binom_inferred_repeats = []
if len(cleaned_trans_associations)!=0:
    for potential_repeat in cleaned_trans_associations:
        if potential_repeat!=primary_junctions_trans:
            # check if the out output constitutes a repeat, if not, discard it so that singleton is only option
            if len(set(extract_ele_list_pairs(potential_repeat))) == 2 and len(potential_repeat) == 2:
                binom_inferred_repeats.append(list(set(extract_ele_list_pairs(potential_repeat))))

            elif len(set(extract_ele_list_pairs(potential_repeat))) == 1 and len(potential_repeat) == 1:
                binom_inferred_repeats.append(list(set(extract_ele_list_pairs(potential_repeat))))

            elif 2* len(set(extract_ele_list_pairs(potential_repeat))) == len(potential_repeat):
                binom_inferred_repeats.append(list(set(extract_ele_list_pairs(potential_repeat))))

#extract junction reference boundaries and store in dictionary
junction_ref_boundaries_dict = return_junction_ref_boundaries(true_junctions)
junction_type_dict = return_junction_inversion_status(true_junctions)

#print junction_ref_boundaries_dict, "here"

WT_content = 0
Petite_content = 0
confident_Petite_content = 0
confident_WT_content = 0
for read_alphabet, read_orientation, read_low_high, read_endpoints, read_juncformat, read_lengths in zip(wholeref_read_alphabet, wholeref_read_alphabet_orientation, wholeref_read_alphabet_num_order, wholeref_read_endpoints, returned_read_juncformat_wholeref, returned_read_lengths_wholeref):

    if len(read_alphabet) ==0:
        WT_content+= read_lengths
        #print read_alphabet, read_endpoints, read_lengths
    #case of wt read but containing an artifact invdup
    elif len(read_alphabet) == 1 and '-' in read_alphabet:
        WT_content+=read_lengths
    #everything else goes to petite structure
    else:
        Petite_content+=read_lengths
        #print read_alphabet, read_endpoints, read_lengths

    #do a version where we only believe junctions are true if they aren't accompanied by an invdup artifact
    if len(read_alphabet) ==0:
        confident_WT_content+= read_lengths
        #print read_alphabet, read_endpoints, read_lengths
    #case of wt read but containing an artifact invdup
    elif '-' in read_alphabet or '*' in read_alphabet:
        confident_WT_content+=read_lengths
    #everything else goes to petite structure
    else:
        confident_Petite_content+=read_lengths
        #print read_alphabet, read_endpoints, read_lengths

    print read_alphabet, read_endpoints, read_lengths

petite_like_frac =  float(Petite_content)/(Petite_content + WT_content)
confident_petite_like_frac = float(confident_Petite_content)/(confident_Petite_content + confident_WT_content)

currDir = os.path.dirname(os.path.realpath(__file__))
folder_name = '/petite_fraction/'
path = currDir + folder_name

if not os.path.exists(path):
    os.mkdir(path)

with open(path + "YPD_Petite_like_fraction.txt", 'a') as f:
    f.write(output_name + "\t" + str(petite_like_frac) + "\t" + str(confident_petite_like_frac) + "\n")


    # if len(read_alphabet)>0:
    #     read_alphabet_nums = [abs(ele) for ele in read_alphabet if isinstance(ele, int)]
    #     read_alphabet_legal_repeats = [ele for ele in read_alphabet_nums if ele in true_junctions_in_repeats]
    #
    #     print read_alphabet
    #     if len(read_alphabet_legal_repeats)>0:
    #
    #         print read_alphabet, read_low_high, read_orientation
    #
    #         query_alignments = get_alignments_from_read(read_alphabet, read_endpoints, true_junctions, read_juncformat)
    #
    #         intersections= [] #at each position, junctions that intersect repeat (primary or alternate, first is primary)
    #         repeat_size = [] #repeat size, which is length of lcs_keys[0] for alternates, len(set(from_primary)) for primary
    #         repeat_junctions = [] #junctions included in detected repeat
    #         estimated_repeat_num = [] #simply minimum count/repeat size for each class of repeats
    #
    #         fractional_nonoverlaps = [] #fraction of nonoverlapping alignments to total ref alignments
    #         fractional_nonoverlaps.append(return_fractional_nonoverlap(query_alignments, primary_alignments))
    #
    #         #populate intersection ,repeat_size list for primary
    #         intersections.append(list(set(primary_aln_key_count[0]) & set(read_alphabet_legal_repeats)))
#
#             if primary_mixed_repeat:
#                 repeat_size.append(len(set(primary_aln_key_count[0])))
#             else:
#                 repeat_size.append(len(primary_aln_key_count[2][0]))
#
#             #populate intersection, repeat_size lists for alternates
#             if len(alternate_aln_key_count)>0:
#                 for alternate, alternate_aln in zip(alternate_aln_key_count, alternate_alignments):
#                     intersections.append(list(set(alternate[0]) & set(read_alphabet_legal_repeats)))
#                     fractional_nonoverlaps.append(return_fractional_nonoverlap(query_alignments, alternate_aln))
#                     repeat_size.append(len(alternate[2][0]))
#
#             #for each entry in intersections, compute estimated_repeat_num present for each possible class
#             for i, inter in enumerate(intersections):
#                 counts = []
#                 for junc in inter:
#                     counts.append(read_alphabet_legal_repeats.count(junc))
#                 if len(counts)!=0:
#                     min_count = min(counts)
#                 else:
#                     min_count=0
#
#                 estimated_repeat_num.append(float(min_count)*len(inter)/repeat_size[i])
#
#             best_matches = np.where(np.array(estimated_repeat_num) == np.array(estimated_repeat_num).max())[0]
#
#             #if there are multiple best matches, take one with lowest nonoverlap fraction
#             if len(best_matches)!=1:
#                 potential_indicies=[]
#                 potential_nonfrac_overlaps = []
#                 for index in best_matches:
#                     potential_indicies.append(index)
#                     potential_nonfrac_overlaps.append(fractional_nonoverlaps[index])
#
#                 #now find corresponding minimum in fractional nonoverlaps, and record
#                 lowest_fracs = np.where(np.array(potential_nonfrac_overlaps) == np.array(potential_nonfrac_overlaps).min())[0]
#                 suggested_indicies = [potential_indicies[i] for i in lowest_fracs]
#
#                 if len(suggested_indicies)==1:
#                     assignment = suggested_indicies[0]
#                 elif len(suggested_indicies) >1 and 0 in suggested_indicies:
#                     assignment = 0
#
#             else:
#                 assignment = best_matches[0]
#
#             #for a special case, if assignment is for primary, but alignment is not consistent with primary,
#             if assignment==0 and fractional_nonoverlaps[0]>0.05:
#
#                 assignment = None
#
#                 # construct inferred primary option here... (if single index make high low etc)
#                 if len(primary_aln_key_count[0])==1 and inferred_present==False:
#                     inferred_present=True
#                     potential_alignment = [junction_ref_boundaries_dict[primary_aln_key_count[0][0]][0][2],
#                                            junction_ref_boundaries_dict[primary_aln_key_count[0][0]][1][2],
#                                            junction_ref_boundaries_dict[primary_aln_key_count[0][0]][2][2],
#                                            junction_ref_boundaries_dict[primary_aln_key_count[0][0]][3][2]]
#
#                     inferred_primary = [primary_aln_key_count[0], potential_alignment,[(str(primary_aln_key_count[0][0]) + 'LH',)],1,0]
#                 elif len(primary_aln_key_count[0])==1:
#                     inferred_primary[-2]+=1
#                     inferred_primary[-1] += 1
#
#             #do incrementing on existing lists
#             if isinstance(assignment, int):
#                 if assignment==0:
#                     primary_aln_key_count[-2] += 1
#                     if '-' in read_alphabet:
#                         primary_aln_key_count[-1] += estimated_repeat_num[assignment]/2.0
#                     else:
#                         primary_aln_key_count[-1] += estimated_repeat_num[assignment]
#                 else:
#                     alternate_aln_key_count[assignment - 1][-2] += 1
#                     if '-' in read_alphabet:
#                         alternate_aln_key_count[assignment - 1][-1] += estimated_repeat_num[assignment]/2.0
#                     else:
#                         alternate_aln_key_count[assignment - 1][-1] += estimated_repeat_num[assignment]
#
#             #if we have a mixed repeat append an example (if '-' not in read, and if estimated_repeat_num > 2)
#             if provide_example_mixed_read and assignment==0:
#                 if len(read_alphabet_legal_repeats)>2:
#                     #print read_alphabet, read_low_high
#                     if (assignment==0) and '-' not in read_alphabet and '*' not in read_alphabet:
#                         #if has legal mixed junctions
#                         #print read_alphabet, read_low_high
#                         if sorted(set([abs(ele) for ele in read_alphabet])) == sorted(set(list(primary_trans_dict.keys()[0]))):
#
#                             #if all junction transitions are legal too
#                             all_transitions_legal=True
#                             for z in range(0, len(read_alphabet)-1):
#                                 transition = [abs(ele) for ele in read_alphabet[z:z+2]]
#                                 if sorted(set(transition))!=sorted(set(list(primary_trans_dict.keys()[0]))):
#                                     all_transitions_legal=False
#
#                             if all_transitions_legal:
#                                 example_mixed_read[tuple(read_low_high)] = [read_orientation, read_alphabet]
#
#
# #print "example test here", example_mixed_read
# if not inferred_present:
#     inferred_primary = {}
#
# #if we have already capture inferred primary in primary list
# if inferred_present:
#     inferred_alignment_start = inferred_primary[1][0]
#     inferred_alignment_end = inferred_primary[1][1]
#
#     for alignment in primary_aln_key_count[1]:
#         if inferred_alignment_start in alignment and inferred_alignment_end in alignment:
#             inferred_primary = {}
#
# if not primary_mixed_repeat:
#     mixed_alignment_count_dict = {}
#     example_mixed_read = {}
# else:
#     example_mixed_read = add_return_ref_locations_dict_example(junction_ref_boundaries_dict, example_mixed_read)
#
# #now iterate through example_mixed_read keys and sort by largest
# if primary_mixed_repeat and len(example_mixed_read.keys())!=0:
#     keys = []
#     key_size = []
#     for key in example_mixed_read.keys():
#         keys.append(key)
#         key_size.append(len(key))
#
#     #now sort, remove keys we don't care about (only keep top 5 longest)
#     sorted_keys, sorted_key_size = (list(t) for t in  zip(*sorted(
#             zip(keys, key_size), key=lambda pair: pair[1], reverse=True )))
#
#     example_mixed_to_remove = sorted_keys[5:]
#     for remove in example_mixed_to_remove:
#         example_mixed_read.pop(remove)
#
# #for a final case, if the span of primary alignment is less than the alternate alignments, make this alternate primary
# # handles 6I2 case... maybe if not robust turn off for rest
# if output_name.split('_')[0]=='6I2':
#     primary_min_refposition = min([min(ele[:2]) for ele in primary_aln_key_count[1]])
#     primary_max_refposition = max([max(ele[:2]) for ele in primary_aln_key_count[1]])
#
#     print primary_min_refposition, primary_max_refposition
#
#     alternates_to_remove = []
#     for alternate in alternate_aln_key_count:
#         alternate_min_refposition = min([min(ele[:2]) for ele in alternate[1]])
#         alternate_max_refposition = max([max(ele[:2]) for ele in alternate[1]])
#
#         if alternate_max_refposition- alternate_min_refposition > primary_max_refposition - primary_min_refposition:
#             primary_aln_key_count = alternate
#             mixed_alignment_count_dict = {}
#             alternates_to_remove.append(alternate)
#
#     for ele in alternates_to_remove:
#         alternate_aln_key_count.remove(ele)
#
#     print primary_aln_key_count
#
# '''STAGING FOR OUTPUT'''
#
# print "__________________________________________________________________"
#
# print "Mixed alignment count dict", mixed_alignment_count_dict, "\n"
#
# print example_mixed_read, "\n"
#
# print "primary transition sign change counts", primary_trans_dict
#
# print inferred_primary, "inferred primary"
#
# print "primary_lcs_suggestions", primary_lcs_suggestion_dict_ref, "\n"
#
# print "primary with counts", primary_aln_key_count, "\n"
#
# print "alternate with counts", alternate_aln_key_count, "\n"
#
# print junction_ref_boundaries_dict, "\n"
#
# print junction_type_dict, "\n"
#
# print junction_counts_dict, "\n"
#
# print "__________________________________________________________________"

#packaging into outputlist

# mixed_alignment_count_dict = {}
# example_mixed_read = {}
# primary_trans_dict = {}
#
# inferred_primary = []
#
# primary_lcs_suggestion_dict_ref = {}
#
# primary_aln_key_count = []
# alternate_aln_key_count = []
#
# output_list = [mixed_alignment_count_dict, example_mixed_read,  primary_trans_dict, inferred_primary,
#                primary_lcs_suggestion_dict_ref, primary_aln_key_count, alternate_aln_key_count,
#                junction_ref_boundaries_dict, junction_type_dict, junction_counts_dict]
#





