import pickle
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec
import numpy as np
import random
import sys
from matplotlib import colors
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch

input_name = "8I1_CP321_30_300_3_1000.pikl"
output_name = input_name[0:-5]

print input_name, "test", output_name
#output_name = sys.argv[1].strip().split('/')[-1][0:-4].split('_')[0] + "_" + sys.argv[1].strip().split('/')[-1][0:-4].split('_')[1]
#print output_name

def get_random_color(pastel_factor = 0.5):
    return [(x+pastel_factor)/(1.0+pastel_factor) for x in [random.uniform(0,1.0) for i in [1,2,3]]]

def color_distance(c1,c2):
    return sum([abs(x[0]-x[1]) for x in zip(c1,c2)])

def generate_new_color(existing_colors,pastel_factor = 0.5):
    max_distance = None
    best_color = None
    for i in range(0,100):
        color = get_random_color(pastel_factor = pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color,c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color

def rotate(l, n):
    return l[-n:] + l[:-n]

currDir = os.path.dirname(os.path.realpath(__file__))
print currDir
#currDir = os.getcwd()
targetDir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
folder_name = '/pikl_pipeline_output_petites/'
path = targetDir + folder_name
print path

del_ins_threshold = input_name.split('.')[0].split('_')[2] #read space vs ref space mismatch required to call junction
min_alignment_length = input_name.split('.')[0].split('_')[3] #mininmum number of bp for alignment to have (removed before junction calling)
min_read_threshold = input_name.split('.')[0].split('_')[4] #minimum number of independent reads that must support one cluster
cluster_merging_eps = input_name.split('.')[0].split('_')[5] #maximum distance between points to cluster in junction space (note invdup is automatic eps)

#load file
per_sample_pikl = []
with open(path + input_name , 'rb') as f:
    per_sample_pikl = pickle.load(f)

'''Labeling inputs and formats for quick referencing'''
mixed_alignment_count_dict = per_sample_pikl[0] # key: LH format, value: [[aln start, aln end, aln start std, aln end std], count]
example_mixed_read_dict = per_sample_pikl[1] # key: LH format, value: [[[LH representation]], [[alphabet representation]], [[aln1, aln end 1, aln1 std, aln end 1 std], [aln2, aln end2, aln2 std, aln 2 end std]]]
primary_trans_dict = per_sample_pikl[2] # key: (transition tuple), value: [count no sign transition, count sign transition]

inferred_primary= per_sample_pikl[3] #[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]
primary_lcs_suggestion_dict = per_sample_pikl[4] # key: (lH repeat format), value : [[orientation rep], [alphabet rep], [[aln1, aln1 end, aln1 std, aln1 end std], [aln2 , aln2 end, aln2 std, aln2 end std]]]
primary_aln_key_count = per_sample_pikl[5] # [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]

alternate_aln_key_count = per_sample_pikl[6] # [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]

junction_ref_boundaries_dict = per_sample_pikl[7] # key: abs(junction num), value: [[LL, LH, L mean], [HL, HH, H mean], [LL std, LH std, L mean std], [HL std, HH std, H mean std]
junction_type_dict = per_sample_pikl[8] # key: abs(junction num), value: True or False, (True meaning inverted junction type, False, noninverted)
junction_counts_dict = per_sample_pikl[9] # key : abs(junction num), value: count (this is for all true junctions whether or not they show up in repeats)

'''Plotting alignment in blocks'''
size = 36

#first generate (10) numbers of distinguishable colours
mycolours = []
names = []
x = 20
for i in range(1,x):
    mycolour = generate_new_color(mycolours,pastel_factor = 0.9)
    mycolours.append(mycolour)


print mycolours[0:4]
#mycolours = ['tab:blue', 'tab:orange', 'tab:pink', 'tab:green']
mycolours = [[0.49599524784617577, 0.998796615225358, 0.9461204803523624], [0.9558576704250413, 0.484576188010616, 0.5291917715595449], [0.6961155955487853, 0.9261837972082733, 0.5082640108662956], [0.8379572484739874, 0.6560311980175444, 0.9630789108265082]]

'''construct primary alignments to colours dict'''
alignment_colour_dict_primary = {}
colour_index = 0

if len(mixed_alignment_count_dict.keys())!=0:
    for key in mixed_alignment_count_dict.keys():
        alignment = mixed_alignment_count_dict[key][0]
        rotated_alignment = rotate(alignment[:2],1)
        rotated_alignment.extend(rotate(alignment[2:],1))

        if rotated_alignment[0]<rotated_alignment[1]:
            for_key = rotated_alignment
        else:
            for_key = alignment

        if tuple(alignment) not in alignment_colour_dict_primary.keys() and tuple(rotated_alignment) not in alignment_colour_dict_primary.keys():
            alignment_colour_dict_primary[tuple(for_key)] = mycolours[colour_index]
            colour_index+=1

else:
    if len(inferred_primary)!=0:

        alignment = inferred_primary[1]
        rotated_alignment = rotate(alignment[:2], 1)
        rotated_alignment.extend(rotate(alignment[2:], 1))

        if rotated_alignment[0] < rotated_alignment[1]:
            for_key = rotated_alignment
        else:
            for_key = alignment

        if tuple(alignment) not in alignment_colour_dict_primary.keys() and tuple(rotated_alignment) not in alignment_colour_dict_primary.keys():
            alignment_colour_dict_primary[tuple(for_key)] = mycolours[colour_index]
            colour_index += 1

        #actually append primary to alternates then..
        alternate_aln_key_count.append(primary_aln_key_count)

    else:
        for alignment in primary_aln_key_count[1]:
            rotated_alignment = rotate(alignment[:2], 1)
            rotated_alignment.extend(rotate(alignment[2:], 1))

            if rotated_alignment[0] < rotated_alignment[1]:
                for_key = rotated_alignment
            else:
                for_key = alignment

            if tuple(alignment) not in alignment_colour_dict_primary.keys() and tuple(rotated_alignment) not in alignment_colour_dict_primary.keys():
                alignment_colour_dict_primary[tuple(for_key)] = mycolours[colour_index]
                colour_index += 1

'''construct alternate alignments to colours dict'''
alignment_colour_dict_alternate = {}

for ele in alternate_aln_key_count:
    print ele
    if not ele[-1] ==0:
        for alignment in ele[1]:

            rotated_alignment = rotate(alignment[:2], 1)
            rotated_alignment.extend(rotate(alignment[2:], 1))

            if rotated_alignment[0] < rotated_alignment[1]:
                for_key = rotated_alignment
            else:
                for_key = alignment

            if tuple(alignment) not in alignment_colour_dict_alternate.keys() and tuple(rotated_alignment) not in alignment_colour_dict_alternate.keys():
                alignment_colour_dict_alternate[tuple(for_key)] = mycolours[colour_index]
                colour_index += 1

#compute number of lines necessary to scale both axes to make arrows the same size
num_lines_first_axis = 2* (len(alignment_colour_dict_primary.keys()) + len(alignment_colour_dict_alternate.keys()))

if len(mixed_alignment_count_dict.keys())!=0:
    num_lines_second_axis= 2 * (len(alternate_aln_key_count) + len(example_mixed_read_dict.keys()) + len(primary_lcs_suggestion_dict.keys()))

else:
    num_lines_second_axis = 2 * (len(alternate_aln_key_count) + 1)


gs = gridspec.GridSpec(2,1,height_ratios=[1,  num_lines_second_axis/float(num_lines_first_axis)])

'''plot all unique alignments on their own axis'''
num_lines = 2* (len(alignment_colour_dict_primary.keys()) + len(alignment_colour_dict_alternate.keys())) + 1
ax2  = plt.subplot(gs[0])
ax2.set_xticklabels([0, 10, 20, 30, 40, 50, 60, 70, 80])
ax2.set_xlim([0,85800])
ax2.set_ylim(0, num_lines)

ax2.set_xlabel('reference mt genome location (kbp)', fontsize=size)

ax2.set_frame_on(False)
ax2.axes.get_yaxis().set_visible(False)
xmin, xmax = ax2.get_xaxis().get_view_interval()
ymin, ymax = ax2.get_yaxis().get_view_interval()
ax2.tick_params(axis='x', labelsize=size)
#ax2.set_title("Unique alignments", fontsize=1.5*size)
ax2.add_line(Line2D((xmin, xmax), (ymin*0.5, ymin*0.5), color='black', linewidth=4))

#handling primary:
last_vertical_index = num_lines

start_primary = last_vertical_index

for aln in alignment_colour_dict_primary.keys():

    alignment = list(aln)
    rotated_alignment = rotate(alignment[:2], 1)
    rotated_alignment.extend(rotate(alignment[2:], 1))

    start = alignment[0]
    end = alignment[1]

    dx = end - start

    body_w = 0.6
    head_w = 1.5 * body_w
    head_l = 0.9 * dx

    colour = alignment_colour_dict_primary[aln]

    if len(mixed_alignment_count_dict.keys())!=0:
        for key in mixed_alignment_count_dict.keys():
            if alignment in mixed_alignment_count_dict[key] or rotated_alignment in mixed_alignment_count_dict[key]:
                mixed_aln_count = mixed_alignment_count_dict[key][-1]

        if abs(dx) > 1000:
            ax2.arrow(start, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True, head_length=1000,
                      head_width=head_w, facecolor=colour, edgecolor='k', linewidth=2)
        else:
            ax2.arrow(start, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                      head_length=head_l,
                      head_width=head_w, facecolor=colour, edgecolor='k', linewidth=2)
        #ax2.text(start + dx + (0.05)*(start + dx), last_vertical_index - 1, "count: %i" %(mixed_aln_count), fontsize=size)
    else:

        if abs(dx)>1000:
            ax2.arrow(start, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True, head_length=1000,
                      head_width=head_w, facecolor=colour, edgecolor='k', linewidth=2)
        else:
            ax2.arrow(start, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True, head_length=head_l,
                      head_width=head_w, facecolor=colour, edgecolor='k', linewidth=2)

    last_vertical_index-=2
end_primary = last_vertical_index
rect = patches.Rectangle((0, end_primary), 85800, abs(end_primary-start_primary), linewidth=2, edgecolor='k', facecolor='none', capstyle='round', joinstyle='round')

#ax2.add_patch(rect)

# if len(mixed_alignment_count_dict.keys())!=0:
#     ax2.text(200, start_primary-1, "Alignments and counts involved\nin primary (mixed) structure", size = size)
# else:
#     ax2.text(200, start_primary - 1, "Alignments involved in primary structure", size=size)

#handling alternate
start_alternate = last_vertical_index
for aln in alignment_colour_dict_alternate.keys():
    alignment = list(aln)

    start = alignment[0]
    end = alignment[1]

    dx = end - start

    body_w = 0.6
    head_w = 1.5 * body_w
    head_l = 0.9 * dx

    colour = alignment_colour_dict_alternate[aln]
    print last_vertical_index
    if abs(dx)>1000:
        ax2.arrow(start, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True, head_length=1000,
                  head_width=head_w, facecolor=colour, edgecolor='k', linewidth=2)
    else:
        ax2.arrow(start, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True, head_length=head_l,
                  head_width=head_w, facecolor=colour, edgecolor='k', linewidth=2)

    last_vertical_index-=2

end_alternate = last_vertical_index

rect = patches.Rectangle((0, end_alternate -2), 85800, abs(end_alternate-start_alternate)+2, linewidth=2, edgecolor='k', facecolor='none')

#ax2.add_patch(rect)
#ax2.text(200, start_alternate -1 , "Alternate repeat alignments", size = size)


#replace patches with rounded versions
new_patches = []
for patch in reversed(ax2.patches):
    bb = patch.get_bbox()
    color=patch.get_facecolor()
    p_bbox = FancyBboxPatch((bb.xmin, bb.ymin),
                        abs(bb.width), abs(bb.height),
                        boxstyle="round,pad=-0.05",
                            ec='k', facecolor='none',
                            mutation_scale=0.00001
                        )
    patch.remove()
    new_patches.append(p_bbox)
for patch in new_patches:
    ax2.add_patch(patch)




'''now plot repeats (irrespective of axis) on next axis'''
ax3  = plt.subplot(gs[1])
ax3.set_xlim([0,85800])

ax3.set_frame_on(False)
ax3.axes.get_yaxis().set_visible(False)
ax3.axes.get_xaxis().set_visible(False)

#ax3.set_title("Structures observed (arbitrary x-axis)", fontsize=1.5*size)



#include possibility for both LCS example repeats and example mixed repeat
if len(mixed_alignment_count_dict.keys())!=0:
    num_lines = 2 * (len(alternate_aln_key_count) + len(example_mixed_read_dict.keys()) + len(primary_lcs_suggestion_dict.keys())) + 3
    ax3.set_ylim(0, num_lines)

    last_vertical_index = num_lines
    start_example_mixed_repeat = last_vertical_index

    last_vertical_index-=0.5

    for key in example_mixed_read_dict.keys():
        # summing span strictly to center example mixed read
        sum = 0
        for aln in example_mixed_read_dict[key][2]:
            sum += abs(aln[0] - aln[1])

        curr_position = (85779 - sum) / 2.0

        junctions_flank_start = example_mixed_read_dict[key][1][0]
        junctions_internal = example_mixed_read_dict[key][1][1:]
        alignments = example_mixed_read_dict[key][2]

        #add text for first flanking junction
        #ax3.text(curr_position, last_vertical_index-0.5, junctions_flank_start, fontsize=size, horizontalalignment='right')

        for i, aln in enumerate(alignments):

            start = aln[0]
            end = aln[1]

            for key in alignment_colour_dict_primary.keys():
                if start in key and end in key:
                    colour = alignment_colour_dict_primary[key]

            dx = end - start
            body_w = 0.6
            head_w = 1.5 * body_w
            head_l = 0.9 * abs(dx)

            if dx < 0:
                curr_position += abs(dx)
                if abs(dx)>1000:
                    ax3.arrow(curr_position, last_vertical_index-1, dx, 0, width=body_w, length_includes_head=True,
                              head_length=1000,
                              head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                else:
                    ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                              head_length=head_l,
                              head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)

            else:
                if abs(dx)>1000:
                    ax3.arrow(curr_position, last_vertical_index-1, dx, 0, width=body_w, length_includes_head=True,
                              head_length=1000,
                              head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                else:
                    ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                              head_length=head_l,
                              head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                curr_position += dx
            if i == len(alignments)-1:
                print "test"
                #ax3.text(curr_position, last_vertical_index-0.5, junctions_internal[i], fontsize=size, horizontalalignment='left')
            else:
                print "test"
                # ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
                #          horizontalalignment='center')
        last_vertical_index-=2



    read_count = primary_aln_key_count[-2]
    monomer_count = primary_aln_key_count[-1]
    if len(example_mixed_read_dict.keys()) != 0:
        # ax3.text(curr_position + 0.1 * curr_position, last_vertical_index + 0.5,
        #          "RC: %i\nMC: %i" % (read_count, monomer_count), fontsize=size)
        end_example_mixed_repeat = last_vertical_index
        rect = patches.Rectangle((0, end_example_mixed_repeat), 85800, abs(end_example_mixed_repeat - start_example_mixed_repeat), linewidth=2,
                             edgecolor='k', facecolor='none')
        #ax3.add_patch(rect)
        #ax3.text(200, start_example_mixed_repeat - 1.2, "Example (mixed)\nprimary structure reads", size=size)

    # print primary lcs suggestion example if it exists
    if len(primary_lcs_suggestion_dict.keys()) != 0:
        start_lcs_primary = last_vertical_index
        last_vertical_index-=0.5


        for key in primary_lcs_suggestion_dict.keys():

            sum = 0

            for aln in primary_lcs_suggestion_dict[key][2]:
                sum += abs(aln[0] - aln[1])

            curr_position = (85779 - sum) / 2.0

            alignments = primary_lcs_suggestion_dict[key][2]
            junctions = primary_lcs_suggestion_dict[key][1]

            # by default, make first alignment forward facing
            start_index_forward = 0
            have_start_index = False
            for i, aln in enumerate(alignments):
                if aln[1] - aln[0] > 0 and not have_start_index:
                    start_index_forward = i
                    have_start_index = True

            alignments = rotate(alignments, start_index_forward)
            rotated_junctions = rotate(junctions, start_index_forward)
            print rotated_junctions, alignments
            junctions_flank_start = rotated_junctions[0]
            junctions_internal = rotated_junctions[1:]
            junctions_internal.append(junctions_flank_start)

            print len(alignments), len(junctions_internal)
            # add text for first flanking junction
            # ax3.text(curr_position, last_vertical_index - 0.5, junctions_flank_start, fontsize=size,
            #          horizontalalignment='right')

            for i, aln in enumerate(alignments):
                if have_start_index:
                    start = aln[0]
                    end = aln[1]
                else:
                    start = min(aln[0], aln[1])
                    end = max(aln[0], aln[1])

                for key in alignment_colour_dict_primary.keys():
                    if start in key and end in key:
                        colour = alignment_colour_dict_primary[key]

                dx = end - start
                body_w = 0.6
                head_w = 1.5 * body_w
                head_l = 0.9 * abs(dx)

                if dx < 0:
                    curr_position += abs(dx)
                    if abs(dx)>1000:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=1000,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    else:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w,
                                  length_includes_head=True,
                                  head_length=head_l,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)


                else:
                    if abs(dx)>1000:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=1000,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    else:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w,
                                  length_includes_head=True,
                                  head_length=head_l,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)

                    curr_position += dx
                if i == len(alignments)-1:
                    print "test"
                    # ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
                    #      horizontalalignment='left')
                else:
                    print "test"
                    # ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
                    #          horizontalalignment='center')

            last_vertical_index-=2
        end_lcs_primary = last_vertical_index

        rect = patches.Rectangle((0, end_lcs_primary), 85800,
                                 abs(end_lcs_primary- start_lcs_primary), linewidth=2,
                                 edgecolor='k', facecolor='none')
        #ax3.add_patch(rect)
        #ax3.text(200, start_lcs_primary - 1.2, "Example detected primary\nstructural repeat unit", size=size)

    #next do alternates
    start_alternate = last_vertical_index
    last_vertical_index-=0.5
    for ele in alternate_aln_key_count:
        if not ele[-1]==0:
            alignments=ele[1]
            junctions = list(ele[2][0])

            junctions = [int(x[:-2]) for x in junctions]

            start_index_forward = 0
            have_start_index = False
            for i, aln in enumerate(alignments):
                if aln[1] - aln[0] > 0 and not have_start_index:
                    start_index_forward = i
                    have_start_index = True

            alignments = rotate(alignments, start_index_forward)
            rotated_junctions = rotate(junctions, start_index_forward)

            junctions_flank_start = rotated_junctions[0]
            junctions_internal = rotated_junctions[1:]
            junctions_internal.append(junctions_flank_start)

            print junctions_internal, alignments
            sum = 0
            for aln in alignments:
                sum += abs(aln[0] - aln[1])
            curr_position = (85779 - sum) / 2.0

            # add text for first flanking junction
            ax3.text(curr_position, last_vertical_index - 0.5, junctions_flank_start, fontsize=size,
                     horizontalalignment='right')

            for i, aln in enumerate(alignments):
                if have_start_index:
                    start = aln[0]
                    end = aln[1]
                else:
                    start = min(aln[0], aln[1])
                    end = max(aln[0], aln[1])

                for key in alignment_colour_dict_alternate.keys():
                    if start in key and end in key:
                        colour = alignment_colour_dict_alternate[key]

                dx = end - start
                body_w = 0.6
                head_w = 1.5 * body_w
                head_l = 0.9 * abs(dx)

                if dx < 0:
                    curr_position += abs(dx)
                    if abs(dx) > 1000:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=1000,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    else:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=head_l,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)

                else:
                    if abs(dx) > 1000:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=1000,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    else:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=head_l,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)

                    curr_position += dx
                if i == len(alignments)-1:

                    ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
                             horizontalalignment='left')
                else:
                    ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
                             horizontalalignment='center')

            read_count = ele[-2]
            monomer_count = ele[-1]
            ax3.text(curr_position + 0.05*curr_position, last_vertical_index - 1.5,
                     "RC: %i\nMC: %i" % (read_count, monomer_count), fontsize=size)
            last_vertical_index -= 2

    end_alternate = last_vertical_index

    rect = patches.Rectangle((0, end_alternate), 85800,
                             abs(end_alternate - start_alternate), linewidth=2,
                             edgecolor='k', facecolor='none', joinstyle='round')
    #ax3.add_patch(rect)
    #ax3.text(200, start_alternate - 0.7, "Alternate repeat units", size=size)

else: #do non mixed primary structure version

    #just primary, followed by alternate
    num_lines = 2 * (len(alternate_aln_key_count) + 1) + 1
    ax3.set_ylim(0, num_lines)
    print num_lines, "num_lines"
    last_vertical_index = num_lines
    start_primary = last_vertical_index
    last_vertical_index-=0.5
    # summing span strictly to center example mixed read
    print inferred_primary

    sum = 0
    if len(inferred_primary)!=0:
        alignments = [inferred_primary[1]]
        monomer_count = inferred_primary[-1]
        read_count = inferred_primary[-2]
    else:
        alignments = primary_aln_key_count[1]
        monomer_count = primary_aln_key_count[-1]
        read_count = primary_aln_key_count[-2]

    junctions = list(primary_aln_key_count[2][0])
    junctions = [x[:-2] for x in junctions]

    for aln in alignments:
        sum += abs(aln[0] - aln[1])

    curr_position = (85779 - sum) / 2.0

    start_index_forward = 0
    have_start_index = False
    for i, aln in enumerate(alignments):
        if aln[1] - aln[0] > 0 and not have_start_index:
            start_index_forward = i
            have_start_index = True

    alignments = rotate(alignments, start_index_forward)
    rotated_junctions = rotate(junctions, start_index_forward)

    junctions_flank_start = rotated_junctions[0]
    junctions_internal = rotated_junctions[1:]
    junctions_internal.append(junctions_flank_start)

    # add text for first flanking junction
    # ax3.text(curr_position, last_vertical_index - 0.5, junctions_flank_start, fontsize=size,
    #          horizontalalignment='right')

    for i, aln in enumerate(alignments):
        if have_start_index:
            start = aln[0]
            end = aln[1]
        else:
            start = min(aln[0], aln[1])
            end = max(aln[0], aln[1])

        for key in alignment_colour_dict_primary.keys():
            if start in key and end in key:
                colour = alignment_colour_dict_primary[key]

        dx = end - start
        body_w = 0.6
        head_w = 1.5 * body_w
        head_l = 0.9 * abs(dx)

        if dx < 0:
            curr_position += abs(dx)
            if abs(dx) > 1000:
                ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                          head_length=1000,
                          head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
            else:
                ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                          head_length=head_l,
                          head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)

        else:
            if abs(dx) > 1000:
                ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                          head_length=1000,
                          head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
            else:
                ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                          head_length=head_l,
                          head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
            curr_position += dx
        # ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
        #          horizontalalignment='center')

    ax3.text(curr_position + 0.05*curr_position, last_vertical_index - 1,
             "RC: %i\nMC: %i" % (read_count, monomer_count), fontsize=size, verticalalignment='center')
    last_vertical_index -= 2
    end_primary = last_vertical_index

    rect = patches.Rectangle((0, end_primary), 85800,
                             abs(end_primary - start_primary), linewidth=2,
                             edgecolor='k', facecolor='none', joinstyle='round')
    #ax3.add_patch(rect)
    ax3.text(200, start_primary - 0.95, "Primary repeat unit", size=size)

    #do alternate
    start_alternate = last_vertical_index
    last_vertical_index-=0.5
    for ele in alternate_aln_key_count:
        if not ele[-1]==0:
            print "***", ele
            alignments = ele[1]
            junctions = list(ele[2][0])
            junctions = [int(x[:-2]) for x in junctions]
            start_index_forward = 0
            have_start_index = False
            for i, aln in enumerate(alignments):
                if aln[1] - aln[0] > 0 and not have_start_index:
                    start_index_forward = i
                    have_start_index = True

            alignments = rotate(alignments, start_index_forward)
            rotated_junctions = rotate(junctions, start_index_forward)

            junctions_flank_start = rotated_junctions[0]

            junctions_internal = rotated_junctions[1:]
            junctions_internal.append(junctions_flank_start)

            print len(junctions_internal), len(alignments)
            print junctions_internal, alignments
            sum=0
            for aln in alignments:
                sum += abs(aln[0] - aln[1])
            curr_position = (85779 - sum) / 2.0

            # add text for first flanking junction
            # ax3.text(curr_position, last_vertical_index - 0.5, junctions_flank_start, fontsize=size,
            #          horizontalalignment='right')

            for i, aln in enumerate(alignments):
                if have_start_index:
                    start = aln[0]
                    end = aln[1]
                else:
                    start = min(aln[0], aln[1])
                    end = max(aln[0], aln[1])

                for key in alignment_colour_dict_alternate.keys():
                    if start in key and end in key:
                        colour = alignment_colour_dict_alternate[key]

                dx = end - start
                body_w = 0.6
                head_w = 1.5 * body_w
                head_l = 0.9 * abs(dx)

                if dx < 0:
                    curr_position += abs(dx)
                    if abs(dx) > 1000:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=1000,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    else:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=head_l,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)

                else:
                    if abs(dx) > 1000:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=1000,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    else:
                        ax3.arrow(curr_position, last_vertical_index - 1, dx, 0, width=body_w, length_includes_head=True,
                                  head_length=head_l,
                                  head_width=head_w, edgecolor='k', facecolor=colour, linewidth=2)
                    curr_position += dx
                # ax3.text(curr_position, last_vertical_index - 0.5, junctions_internal[i], fontsize=size,
                #          horizontalalignment='center')

            read_count = ele[-2]
            monomer_count = ele[-1]
            ax3.text(curr_position  + 0.05 * curr_position, last_vertical_index - 1,
                     "RC: %i\nMC: %i" % (read_count, monomer_count), fontsize=size, verticalalignment='center')

            last_vertical_index -= 2
    end_alternate = last_vertical_index

    rect = patches.Rectangle((0, end_alternate), 85800,
                             abs(end_alternate - start_alternate), linewidth=2,
                             edgecolor='k', facecolor='none', joinstyle='round')
    #ax3.add_patch(rect)
    ax3.text(200, start_alternate - 0.95, "Alternate repeat units", size=size)

#replace patches with rounded versions
new_patches = []
for patch in reversed(ax3.patches):
    bb = patch.get_bbox()
    color=patch.get_facecolor()
    p_bbox = FancyBboxPatch((bb.xmin, bb.ymin),
                        abs(bb.width), abs(bb.height),
                        boxstyle="round,pad=-0.04",
                            ec='k', facecolor='none',
                            mutation_scale=0.0001
                        )
    patch.remove()
    new_patches.append(p_bbox)
for patch in new_patches:
    ax3.add_patch(patch)

#mpl.rcParams['axes.linewidth'] = 20 #set the value globally
fig = plt.gcf()
fig.set_size_inches(14, 20)
fig.tight_layout(w_pad=0.4, h_pad=3)
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(10)
    ax3.spines[axis].set_linewidth(0.5)

fig.savefig(output_name + ".png", dpi=300, format='PNG', transparent=True)