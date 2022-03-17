
This repository contains preprocessed data and analysis code from https://doi.org/10.1101/2021.11.15.468706.

Each folder contains the code for analysis and plotting of each figure in paper.  It also includes folders of PIKL files that are the output of the structural inference pipeline described in the paper. 

To perform your own analysis from this preprocessed data, use the PIKL files. Note that not all indicies of the PIKL are nonzero depending on the sample. The PIKL data structure is below, where aln stands for alignment:

index 0 - mixed alignment count dictionary, format:

    key: LH juntion format, value: [[aln start, aln end, aln start std, aln end std], count]

index 1 - example mixed dictionary, format:

    key: LH junction format, value: [[[LH representation]], [[alphabet representation]], [[aln1, aln end 1, aln1 std, aln end 1 std], [aln2, aln end2, aln2 std, aln 2 end std]]]

index 2 - primary junction transition dictionary, format:

    key: (transition tuple), value: [count no sign transition, count sign transition]

index 3 - inferred primary structure list if nonzero, format:

    [[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]

index 4 - primary repeat suggested from longest common prefix detection, format:

    key: (lH repeat format), value : [[orientation rep], [alphabet rep], [[aln1, aln1 end, aln1 std, aln1 end std], [aln2 , aln2 end, aln2 std, aln2 end std]]]

index 5 - primary alignments and counts list, format:

    [[set of abs junctions involved], [[[aln1, aln1 end, aln1 std, aln1 end std], [aln2, aln2 end, aln2 std, aln2 end std]], [(lH key)], molecular count, monomer count]

index 6 - alternate alignments and counts list, format:

    [[[set of abs junctions involved], [[aln1, aln1 end, aln1 std, aln1 end std], [(lH key)], molecular count, monomer count]]

index 7 - dictionary of all junctions in reference coordinates, format:

    key: abs(junction num), value: [[LL, LH, L mean], [HL, HH, H mean], [LL std, LH std, L mean std], [HL std, HH std, H mean std]

index 8 - dictionary of all junctions types corresponding to 7, format:

    key: abs(junction num), value: True or False, (True meaning inverted junction type, False, noninverted)

index 9 - dictionary of all junctions counts corresponding to 7&8, format:

    key : abs(junction num), value: count (this is for all true junctions whether or not they show up in repeats)







