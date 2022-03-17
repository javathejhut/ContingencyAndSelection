import pysam
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import math
from matplotlib.lines import Line2D
import glob
import os

'''Must change for nuclear vs mt for annotated features filename and header, chr length etc'''
header = 6 # how many lines before data actually starts?
featurefile = "All Annotated Sequence Features-chrmt-1..85779_culled.sqn"

features = {}
raw_lines = []
line_feature_pos = []
bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
print sys.argv[1].strip().split('/')[-1][0:-4]
output_name = sys.argv[1]
complete_read_list = []  # all entries in BAM, included supplementary and secondary reads
unique_read_list = []  # only primary reads in BAM

with open(featurefile, 'r') as annotated:
    # loop to get feature line positions in file
    for pos, line in enumerate(annotated):
        # if we are starting a new feature
        if not line.startswith("\t") and (pos >= header):
            line_feature_pos.append(pos)
        if (pos >= header):  # if we are outside the whole "chromosome" annotation which we definitely don't want
            raw_lines.append(line)

def sequin_format(argument):
    if len(argument) != 0 and len(argument[0].strip().split()) <= 2:
        return argument[0].strip().split()[1]
    else:
        return "N/A"


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

# construct read lists
for entry in bamfile:
    lor = entry.query_length
    if lor > 0:
        if entry.is_unmapped == False and entry.is_supplementary == False and entry.is_secondary == False:
            if not "offset" in entry.query_name.strip().split('_'):
                unique_read_list.append([entry, entry.query_length])
        if entry.is_unmapped == False and entry.is_secondary == False:
            complete_read_list.append([entry, entry.query_length])

total = len(unique_read_list)  # number of unique reads
tail_read_count = total - int(0.90 * total)  # how many unique read names we want in the tail of the size distribution

# list of unique primary reads
sorted_unique_read_list = sorted(unique_read_list, key=lambda x: x[1], reverse=True)
tail_primary_reads = sorted_unique_read_list[:tail_read_count]

# for the tail of these primary reads (which just means unique) get the positions where portions of this read are mapped for primary + supplementary reads
read_loc_dict = {}
alignment_gaps_dict = {}
inverted_junction_signals_dict = {}
normal_junction_signals_dict = {}

for read_p in tail_primary_reads:
    subalignments = []
    read_length = read_p[0].infer_read_length()
    for read in complete_read_list:
        # if a non clipped segment
        if read_p[0].query_name == read[0].query_name and read[0].mapping_quality >= 20:
            clips = []
            clip_positions = []
            cigar_tuples = read[0].cigartuples
            blocks = read[0].get_blocks()
            matches_insertions_wholeread = []
            temp_list = []

            if read[0].is_reverse:
                cigar_tuples = cigar_tuples[::-1]
                blocks = blocks[::-1]

            for i, tup in enumerate(cigar_tuples):
                # if either soft/hard clipped (handles primary read which is only soft-clipped)
                if tup[0] == 4 or tup[0] == 5:
                    clips.append(int(tup[1]))
                    clip_positions.append(i)

            if len(clips) == 2:
                qstart = clips[0]  # starting position of segment in read coordinates
                qend = read_length - clips[-1]
            elif len(clips) == 1:
                if clip_positions[0] == 0:
                    qstart = clips[0]  # starting position of segment in read coordinates
                    qend = read_length
                else:
                    qstart = 0  # starting position of segment in read coordinates
                    qend = read_length - clips[-1]
            else:
                qstart = 0  # starting position of segment in read coordinates
                qend = read_length

            ref_start = read[0].reference_start
            ref_end = read[0].reference_end

            subalignments.append([qstart, qend, ref_start, ref_end, read[0].is_reverse])

            '''construct lists that tell us how many matches there are and the number of Insertions between each match'''
            for index in range(len(cigar_tuples) - 1):

                # if a match, append tuple to list
                if cigar_tuples[index][0] == 0:
                    temp_list.append(cigar_tuples[index])
                # if insertion append tuple to list
                elif cigar_tuples[index][0] == 1:
                    temp_list.append(cigar_tuples[index])
                # if next element matches, append existing matches and insertions to master list
                if cigar_tuples[index + 1][0] == 0 and len(temp_list) != 0:
                    matches_insertions_wholeread.append(temp_list)
                    temp_list = []

                # if we are at last match
                if index + 1 == len(cigar_tuples) - 2:
                    matches_insertions_wholeread.append(temp_list)

            '''now construct the equivalent of blocks but in read coordinates'''
            blocks_in_read = []
            start = qstart
            end = matches_insertions_wholeread[0][0][1] + start

            # first matching portion
            blocks_in_read.append([start, end])

            for ele in matches_insertions_wholeread[1:]:

                # either tuples in list are [match tuple, insertion tuple]
                if len(ele) == 2:
                    num_insertions = ele[1][1]
                    num_matches = ele[0][1]
                    start = end + num_insertions
                    end = start + num_matches

                # or just [match tuple] (likely end of matched read)
                elif len(ele) == 1:
                    num_matches = ele[0][1]
                    start = end
                    end = start + num_matches

                blocks_in_read.append([start, end])

            '''add all this info to a dictionary for each segment of read'''
            if read_p[0].query_name not in read_loc_dict.keys():
                read_loc_dict[read_p[0].query_name] = []


                # format: read,  is_reverse?, q_alignmentstart, q_alignmentend, reference_start, reference_end, blocks, read blocks
                read_loc_dict[read_p[0].query_name].append(
                    [read[0], read[0].is_reverse, qstart, qend, read[0].reference_start,
                     read[0].reference_end,
                     blocks, blocks_in_read])
            else:


                # format: read,  is_reverse?,q_alignmentstart, q_alignmentend, reference_start, reference_end, blocks, read blocks
                read_loc_dict[read_p[0].query_name].append(
                    [read[0], read[0].is_reverse, qstart, qend, read[0].reference_start,
                     read[0].reference_end,
                     blocks, blocks_in_read])

        # if this is from a clipped region, make sure to include offsets
        elif read_p[0].query_name == read[0].query_name.split("offset")[0][:-1] and read[0].mapping_quality >= 20:
            clips = []
            clip_positions = []
            cigar_tuples = read[0].cigartuples
            blocks = read[0].get_blocks()
            matches_insertions_wholeread = []
            temp_list = []

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

                qstart = clips[0] + start_offset  # starting position of segment in read coordinates
                qend = end_offset - clips[-1]

            elif len(clips) == 1:

                if clip_positions[0] == 0:
                    qstart = clips[0] + start_offset  # starting position of segment in read coordinates
                    qend = end_offset

                else:
                    qstart = start_offset  # starting position of segment in read coordinates
                    qend = end_offset - clips[-1]
            else:
                qstart = start_offset  # starting position of segment in read coordinates
                qend = end_offset

            ref_start = read[0].reference_start
            ref_end = read[0].reference_end

            subalignments.append([qstart, qend, ref_start, ref_end, read[0].is_reverse])

            '''construct lists that tell us how many matches there are and the number of Insertions between each match'''
            for index in range(len(cigar_tuples) - 1):

                # if a match, append tuple to list
                if cigar_tuples[index][0] == 0:
                    temp_list.append(cigar_tuples[index])
                #if insertion append tuple to list
                elif cigar_tuples[index][0] == 1:
                    temp_list.append(cigar_tuples[index])
                #if next element matches, append existing matches and insertions to master list
                if cigar_tuples[index + 1][0] == 0 and len(temp_list) != 0:
                    matches_insertions_wholeread.append(temp_list)
                    temp_list = []

                # if we are at last match
                if index + 1 == len(cigar_tuples) - 2:
                    matches_insertions_wholeread.append(temp_list)

            #special case is we have clipped tuple then match, or match then clipped tuple
            blocks_in_read = []
            if len(cigar_tuples)<=2:

                start = qstart
                for tup in cigar_tuples:
                    if tup[0] == 0:
                        #print cigar_tuples, blocks
                        blocks_in_read.append([start, start + tup[1]])

            else:
                '''now construct the equivalent of blocks but in read coordinates'''
                start = qstart

                end = matches_insertions_wholeread[0][0][1] + start

                #first matching portion
                blocks_in_read.append([start, end])

                for ele in matches_insertions_wholeread[1:]:

                    #either tuples in list are [match tuple, insertion tuple]
                    if len(ele) == 2:
                        num_insertions = ele[1][1]
                        num_matches = ele[0][1]
                        start = end + num_insertions
                        end = start + num_matches

                    #or just [match tuple] (likely end of matched read)
                    elif len(ele) == 1:
                        num_matches = ele[0][1]
                        start = end
                        end = start + num_matches

                    blocks_in_read.append([start, end])


            #print len(blocks), len(blocks_in_read)
            '''add all this info to a dictionary for each segment of read'''
            if read_p[0].query_name not in read_loc_dict.keys():
                read_loc_dict[read_p[0].query_name] = []

                #format: read,  is_reverse?, q_alignmentstart, q_alignmentend, reference_start, reference_end, blocks, read blocks
                read_loc_dict[read_p[0].query_name].append([read[0], read[0].is_reverse, qstart, qend, read[0].reference_start, read[0].reference_end,
                                                          blocks, blocks_in_read])
            else:

                # format: read,  is_reverse?,q_alignmentstart, q_alignmentend, reference_start, reference_end, blocks, read blocks
                read_loc_dict[read_p[0].query_name].append([read[0], read[0].is_reverse, qstart, qend, read[0].reference_start, read[0].reference_end,
                                                          blocks, blocks_in_read])

    sorted_subalignments = sorted(subalignments, key=lambda x:x[0])
    if len(sorted_subalignments)!=0:
        alignment_gaps = []
        #for first and last subalignment, compare the start and end respectively to start and end of total read (0 and read_length)
        read_start = 0
        read_end = read_length
        first_subalignment_start = sorted_subalignments[0][0]
        last_subalignment_end = sorted_subalignments[-1][1]

        alignment_gaps.append([read_start, first_subalignment_start])
        alignment_gaps.append([last_subalignment_end, read_end])

        #for all subalignments record missing portions of read between start and end of adjacent alignments
        for i in range(len(sorted_subalignments)-1):
            current_alignment = sorted_subalignments[i]
            next_alignment = sorted_subalignments[i+1]

            current_end = current_alignment[1]
            next_start = next_alignment[0]

            if (next_start - current_end) > 100:
                alignment_gaps.append([current_end, next_start])
        if read_p[0].query_name not in alignment_gaps_dict.keys():

            alignment_gaps_dict[read_p[0].query_name] = []
            alignment_gaps_dict[read_p[0].query_name] = alignment_gaps
        else:
            alignment_gaps_dict[read_p[0].query_name] = alignment_gaps

        #segment sorted_subalignments into normal_junction_signals or inverted_junction_signals
        normal_junction_signals = []
        inverted_junction_signals = []

        for i in range(len(sorted_subalignments)):
            curr_orientation = sorted_subalignments[i][4]

            if curr_orientation==True:
                if read_p[0].query_name not in inverted_junction_signals_dict.keys():
                    inverted_junction_signals_dict[read_p[0].query_name] = []

                inverted_junction_signals_dict[read_p[0].query_name].append(sorted_subalignments[i][0])
                inverted_junction_signals_dict[read_p[0].query_name].append(sorted_subalignments[i][1])

            elif curr_orientation==False:
                if read_p[0].query_name not in normal_junction_signals_dict.keys():
                    normal_junction_signals_dict[read_p[0].query_name] = []
                normal_junction_signals_dict[read_p[0].query_name].append(sorted_subalignments[i][0])
                normal_junction_signals_dict[read_p[0].query_name].append(sorted_subalignments[i][1])

#print read_loc_dict['1fe10929-b8d6-4f4e-b10a-228ac2cab188'][0]

print "-------------------------------------------------------"
read_coding_regions = {}

print len(read_loc_dict.keys())

# now loop through features, and ask what features are within windows where query maps to reference
for read_name in read_loc_dict.keys():

    if len(read_loc_dict[read_name])!=0:
        read_length = max([read[0].infer_read_length() for read in read_loc_dict[read_name]])

        for feature_key in features.keys():
            # replace gene with name for key in features_in_read dict if possible
            feature_gene = features[feature_key][1]
            feature_name = features[feature_key][3]
            feature_start = features[feature_key][5]
            feature_end = features[feature_key][6]

            if feature_gene != 'N/A':
                identifier = feature_gene
            else:
                identifier = feature_name

            #populate list of features in a given read
            if read_name not in read_coding_regions.keys():
                read_coding_regions[read_name] = []

            #for each individual alignment segment (all from one read)
            for segment in read_loc_dict[read_name]:
                blocks_in_ref = segment[6]
                blocks_in_read = segment[7]

                is_reversed = segment[1]

                for index in range(len(blocks_in_ref)):
                    block = blocks_in_ref[index]
                    block_read = blocks_in_read[index]

                    if is_overlapped(block[0], block[1], features[feature_key][5], features[feature_key][6]):

                        q_start = block_read[0]
                        q_end = block_read[1]

                        r_start = block[0]
                        r_end = block[1]

                        #some origins are expected for some reason to be seen in one strand vs other according to SGD
                        #ORI1, ORI8 etc not sure why this is but will have to be careful about forward/reverese strand
                        #mapping because of this
                        feature_loc_in_ref = [int(features[feature_key][5]), int(features[feature_key][6])]

                        f_start = min(feature_loc_in_ref)
                        f_end = max(feature_loc_in_ref)

                        #if left section of feature is cut off
                        if f_start < r_start <= f_end and r_start <= f_end < r_end:
                            q_feature_start = q_start
                            q_feature_end = q_end - (r_end - f_end)

                        #if right section of feature is cut off
                        elif r_start<f_start<=r_end and f_start <=r_end< f_end:
                            q_feature_start = q_end - (r_end-f_start)
                            q_feature_end = q_end

                        #if feature is fully within read:
                        elif r_start<f_start< r_end and r_start<f_end< r_end:
                            q_feature_start = q_start + f_start - r_start
                            q_feature_end = q_end - (r_end - f_end)

                        #if read segment is smaller than feature and contained within it
                        elif f_start<=r_start<f_end and f_start<r_end<=f_end:
                            # print "amihere?"
                            q_feature_start =  q_start
                            q_feature_end = q_end

                        #shouldn't be here, or overlap is only one base
                        else:
                            q_feature_start = q_start
                            q_feature_end = q_end

                        read_coding_regions[read_name].append([identifier, is_reversed, q_feature_start, q_feature_end, int(feature_start), int(feature_end), read_length])

def spiral_plot(readname, path):
    print readname

    final_coding_regions =  read_coding_regions[readname]
    gaps = alignment_gaps_dict[readname]
    if readname in normal_junction_signals_dict.keys():
        normal_junction_signals = normal_junction_signals_dict[readname]
    else:
        normal_junction_signals = []
    if readname in inverted_junction_signals_dict.keys():
        inverted_junction_signals = inverted_junction_signals_dict[readname]
    else:
        inverted_junction_signals = []

    #print final_coding_regions

    read_length = read_coding_regions[readname][0][-1]
    numrepeats = 8 #rounded based on TideHunter output (must divide as an integer with numbases)

    numbases = int(numrepeats * round(read_length / numrepeats))

    revolutions = int(2*numrepeats)
    #numcells = int(2.2*numbases)
    numcells = int(2 * numbases)
    min_r = 1/2.0
    max_r = 21/2.0

    alpha = (revolutions)/((max_r-min_r))

    #plot spiral backbone
    fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(20,20))
    #theta = np.arange(revolutions * 2 * np.pi + np.pi, (revolutions / 2) * 2 * np.pi + np.pi, -0.0001)
    theta = np.arange(revolutions * 2 * np.pi + np.pi, (revolutions / 2) * 2 * np.pi -np.pi, -0.0001)

    r = theta/(alpha*2*np.pi)
    ax.set_axis_off()
    ax.grid(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    '''Constructing legend and colour map'''
    # generate our "reference" which is the order in which we expect to see the genes in sorted_list
    reference_inread = []
    reference = []
    reference_starts = []
    reference_starts_inread = []
    feature_bars = {}
    for key in features.keys():
        feature_gene = features[key][1]
        feature_name = features[key][3]
        feature_start = min([int(features[key][5]), int(features[key][6])])

        for ele in final_coding_regions:
            if (feature_name == ele[0] or feature_gene == ele[0]) and ele[0] not in reference_inread:

                reference_inread.append(ele[0])
                reference_starts_inread.append(feature_start)

        if feature_gene != 'N/A':
            identifier = feature_gene
        else:
            identifier = feature_name

        feature_bars[identifier] = [min(feature_start, feature_end), max(feature_start, feature_end)]

    reference_starts_inread, reference_inread = (list(t) for t in zip(*sorted(zip(reference_starts_inread, reference_inread))))

    #create a reference separate from what is present in read
    feature_names = []
    feature_locations = []
    reference_starts = []

    for key in feature_bars:
        feature_names.append(key)
        feature_locations.append(feature_bars[key])
        reference_starts.append(feature_bars[key][0])

    reference_starts, reference = (list(t) for t in zip(*sorted(zip(reference_starts, feature_names))))

    genes_present = []
    for ele in final_coding_regions:
        if ele[0] not in genes_present:
            genes_present.append(ele[0])

    # print reference
    # print reference_starts
    # print genes_present

    n = len(reference)
    colours = cm.rainbow(np.linspace(0,1,n))

    legend_lines = [Line2D([0],[0], color=colours[reference.index(gene)], lw=10) for gene in reference_inread]

    ax.legend(legend_lines, reference_inread, prop={'size':45}, loc="lower right")
    ax.set_title(str(readname)[0:5] + " read length: " + str(read_length) + " bp", size=45)
    ax.plot(theta, r, linewidth=1, color='k');

    #basepair location labels (start and end)
    delta_theta = 2*np.pi*revolutions/(numcells)
    theta_top = np.arange(np.pi + numcells*(delta_theta), np.pi + (numcells-1)*(delta_theta), -0.0001)
    theta_bottom = np.arange(numbases*(delta_theta) - np.pi, (numbases-1)*(delta_theta) - np.pi, -0.0001)
    r_top = max(theta_top/(alpha*2*np.pi))
    r_bottom = min(theta_bottom/(alpha*2*np.pi))

    plt.text(np.pi, r_top*1.22, "0 bp", fontsize=45)
    plt.text(np.pi, r_bottom*0.95, str(read_length) + " bp", fontsize=45)

    testing_start = []
    testing_end = []
    for ele in final_coding_regions:

        start_position = int(ele[2])
        end_position = int(ele[3])

        start = numcells -  start_position
        end = numcells -  end_position

        testing_start.append(start)
        testing_end.append(end)

        theta_top = np.arange(np.pi + start*(delta_theta), np.pi + end*(delta_theta), -0.0001)
        theta_bottom = np.arange(start*(delta_theta) - np.pi, end*(delta_theta) - np.pi, -0.0001)
        r_top = theta_top/(alpha*2*np.pi)
        r_bottom = theta_bottom/(alpha*2*np.pi)

        if ele[1] == False:  # forward
            ax.fill_between(theta_top, r_bottom, r_top, color=colours[reference.index(ele[0])], alpha=0.5)
        else: #reverse
            ax.fill_between(theta_top, r_bottom, r_top, color=colours[reference.index(ele[0])], hatch='.', alpha=0.5)

    #gaps in alignment, which we will make black to make things clear
    for ele in gaps:

        start_position = int(ele[0])
        end_position = int(ele[1])

        start = numcells -  start_position
        end = numcells -  end_position

        testing_start.append(start)
        testing_end.append(end)

        theta_top = np.arange(np.pi + start*(delta_theta), np.pi + end*(delta_theta), -0.0001)
        theta_bottom = np.arange(start*(delta_theta) - np.pi, end*(delta_theta) - np.pi, -0.0001)
        r_top = theta_top/(alpha*2*np.pi)
        r_bottom = theta_bottom/(alpha*2*np.pi)


        ax.fill_between(theta_top, r_bottom, r_top, color='k', alpha=1.0)

    '''Plotting junctions'''
    for ele in normal_junction_signals:
        start_position = ele - 2
        end_position = ele + 2

        start = numcells - start_position
        end = numcells - end_position

        testing_start.append(start)
        testing_end.append(end)

        theta_top = np.arange(np.pi + start * (delta_theta), np.pi + end * (delta_theta), -0.0001)
        theta_bottom = np.arange(start * (delta_theta) - np.pi, end * (delta_theta) - np.pi, -0.0001)
        r_top = theta_top / (alpha * 2 * np.pi)
        r_bottom = theta_bottom / (alpha * 2 * np.pi)

        ax.fill_between(theta_top, r_bottom, r_top, color='r', alpha=1.0)

    for ele in inverted_junction_signals:
        start_position = ele - 2
        end_position = ele + 2

        start = numcells - start_position
        end = numcells - end_position

        testing_start.append(start)
        testing_end.append(end)

        theta_top = np.arange(np.pi + start * (delta_theta), np.pi + end * (delta_theta), -0.0001)
        theta_bottom = np.arange(start * (delta_theta) - np.pi, end * (delta_theta) - np.pi, -0.0001)
        r_top = theta_top / (alpha * 2 * np.pi)
        r_bottom = theta_bottom / (alpha * 2 * np.pi)

        ax.fill_between(theta_top, r_bottom, r_top, color='b', alpha=1.0)

    plt.gcf()

    plt.savefig(path + str(readname) +'_blocks.png')

#sort keys of dict above by read length so that we can inspect largest reads visually
length_keys = []
length_values = []
for key in read_coding_regions.keys():
    if len(read_coding_regions[key])>0:
        length_keys.append(key)
        length_values.append(int(read_coding_regions[key][0][6]))

length_values, length_keys = (list(t) for t in zip(*sorted(zip(length_values, length_keys))))

currDir = os.path.dirname(os.path.realpath(__file__))
folder_name = '/' + sys.argv[1].strip().split('/')[-1][0:-4] + '/'
path = currDir + folder_name

if not os.path.exists(path):
    os.makedirs(path)

for readname in length_keys[-30:]:

    spiral_plot(readname, path)
