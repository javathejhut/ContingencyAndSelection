import matplotlib.pyplot as plt
import glob
import numpy as np
import pickle
import math
import os

def rotate(l, n):
    return l[-n:] + l[:-n]

mixed_samples = ['3I1', '6I1', '17I1', '23I1', '8I1', '13I1', '19I1']

#start by importing exp fit parameters from mixed sample read length distributions
exp_fit_parameters_dict = {}
with open("readlength_exp_fit_parameters_petites.txt", 'r') as fparameters:
	for line in fparameters:
		name = line.strip().split()[0]
		location = float(line.strip().split()[1])
		scale = float(line.strip().split()[2])

		if name in mixed_samples:
			exp_fit_parameters_dict[name] = [location, scale]

#next import pikls
alignment_dict_alternate = {}
alignment_dict_primary = {}
parent_path = os.path.normpath(os.path.join(os.getcwd(), os.pardir))
print(parent_path + '/pikl_pipeline_output_petites/*.pikl')

for sample_pikl_file in glob.glob(parent_path + '/pikl_pipeline_output_petites/*.pikl'):
	pikl_name= sample_pikl_file.strip().split('/')[-1].split('_')[0]
	print(pikl_name)

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

		'''construct primary alignments '''
		alignment_dict_primary[pikl_name] = []

		if len(mixed_alignment_count_dict.keys()) != 0:

			for key in mixed_alignment_count_dict.keys():
				alignment = mixed_alignment_count_dict[key][0]
				count = mixed_alignment_count_dict[key][1]

				alignment_dict_primary[pikl_name].append([alignment, count])

		else:
			if len(inferred_primary) != 0:

				alignment = inferred_primary[1]
				# actually append primary to alternates then..
				alternate_aln_key_count.append(primary_aln_key_count)

				alignment_dict_primary[pikl_name].append(alignment)
			else:
				for alignment in primary_aln_key_count[1]:
					alignment_dict_primary[pikl_name].append(alignment)


		'''construct alternate alignments '''
		alignment_dict_alternate[pikl_name] = []

		for ele in alternate_aln_key_count:
			for alignment in ele[1]:
				alignment_dict_alternate[pikl_name].append(alignment)

#now parse dictionaries for 1st sample run, comparing median read length and maximum alignment length

mixed_sample_names = []
mixed_alignment_lengths = []
all_samples_aln_1_frac = []
all_samples_aln_2_frac = []
all_samples_aln_3_frac = []
all_samples_aln_4_frac = []

all_samples_1_count = 0
all_samples_2_count = 0
all_samples_3_count = 0
all_samples_4_count = 0

for key_alignment in alignment_dict_primary.keys():
	for key_fit_params in exp_fit_parameters_dict.keys():
		if key_alignment == key_fit_params:
			mu = float(exp_fit_parameters_dict[key_fit_params][0])  #location of exponential
			beta = float(exp_fit_parameters_dict[key_fit_params][1]) #scale of exponential

			print mu, beta
			alignment_lengths = [abs(ele[0][1] - ele[0][0]) for ele in alignment_dict_primary[key_alignment]]
			alignment_counts = [ele[1] for ele in alignment_dict_primary[key_alignment]]

			alignment_lengths, alignment_counts = (list(t) for t in zip(*sorted(zip(alignment_lengths, alignment_counts), reverse=True)))

			# create normalized counts
			normalized_counts = []
			for i, count in enumerate(alignment_counts):
				prob_read_to_contain = 1.0- (alignment_lengths[i]/float(sum(alignment_lengths)))
				normalized_counts.append(count/math.exp(-1 *(alignment_lengths[i] - mu)/beta) * prob_read_to_contain)

			#create probability normalized fractions
			prob_normalized_fractions = []
			for norm_count in normalized_counts:
				prob_normalized_fractions.append(norm_count/sum(normalized_counts))

			#print alignment_counts, normalized_counts, prob_normalized_fractions

			all_samples_aln_1_frac.append(prob_normalized_fractions[0])
			all_samples_aln_2_frac.append(prob_normalized_fractions[1])
			all_samples_aln_3_frac.append(prob_normalized_fractions[2])
			all_samples_aln_4_frac.append(prob_normalized_fractions[3])

			all_samples_1_count+=alignment_counts[0]
			all_samples_2_count += alignment_counts[1]
			all_samples_3_count += alignment_counts[2]
			all_samples_4_count += alignment_counts[3]

to_plot = [np.mean(all_samples_aln_1_frac), np.mean(all_samples_aln_2_frac), np.mean(all_samples_aln_3_frac), np.mean(all_samples_aln_4_frac)]
to_plot_stdev = [np.std(all_samples_aln_1_frac), np.std(all_samples_aln_2_frac), np.std(all_samples_aln_3_frac), np.std(all_samples_aln_4_frac)]

print all_samples_1_count
print all_samples_2_count
print all_samples_3_count
print all_samples_4_count

N = len(to_plot)

ind = np.arange(N)
print ind
print len(ind)
width = 0.5

fig, ax = plt.subplots()

rects = ax.bar(ind, to_plot, width, color=[[0.9558576704250413, 0.484576188010616, 0.5291917715595449],[0.49599524784617577, 0.998796615225358, 0.9461204803523624], [0.6961155955487853, 0.9261837972082733, 0.5082640108662956], [0.8379572484739874, 0.6560311980175444, 0.9630789108265082]], yerr=to_plot_stdev, align='edge')

size = 36

ax.set_ylabel('Average fraction',fontsize=size)
ax.set_xlabel('Alignment length rank', fontsize=size)
ax.set_title('Proportions of alignments across 7 mixed samples\n(corrected for sampling bias)', fontsize=size)
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(['4', '3', '2', '1'], fontsize=size)

#ax.legend((rects[0], rects_chrI[0]), ('chrM', 'chrI'), loc='upper right', fontsize=size)

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.02*height,
                '%2.2f' % (height),
                ha='center', va='bottom', fontsize=size-10)

autolabel(rects)

#plt.legend(loc='upper left')
#plt.legend(loc='upper left', prop={'size': size})
ax.tick_params(axis='both', which='major', labelsize=size)
ax.tick_params(axis='both', which='minor', labelsize=size)

fig = plt.gcf()
#fig.tight_layout()
fig.set_size_inches(16,10)
#plt.show(fig)
#fig = plt.gcf()
#fig.subplots_adjust(bottom=0.5)
fig.savefig("mixed_alignment_proportions.png", DPI=600)






