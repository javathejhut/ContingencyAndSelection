import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys
import time

start_time = time.time()

coverage = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		coverage.append(abs(int(line.split()[2])))
	#print length
	
#for colour map
c = []         #array of colours

#with open("chrM.fa", 'r') as ref:
	
	# #create reference array
	# for line in ref:
	# 	for char in line.strip():
	# 		if char == 'A':
	# 			c.append('g')
	# 		if char == 'T':
	# 			c.append('g')
	# 		else:
	# 			c.append('k')

binwidth = 1000
N = len(coverage)
x = range(N)
#width = 2/3.0
width = 2/3.0
#cmap, norm = plt.colors.from_levels_and_colors([0, 2, 5, 6], ['red', 'green', 'blue'])
plt.bar(x,coverage,width,snap=False, color='black')

print "done bar"

#ax.scatter(x, y, norm=norm)

#plt.hist(coverage, bins=range(0, 85799 + binwidth, binwidth))
size = 30
ax = plt.axes()
ax.margins(0.05)
line_pos = 1.05*max(coverage)
text_pos = 1.01*line_pos
ax.set_xticklabels([0, 10, 20, 30, 40, 50, 60, 70, 80])
plt.xticks([0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000])
plt.yticks([250, 500, 750, 1000, 1250, 1500])
text_pos_lower = 0.96* line_pos

plt.xlabel('reference mt genome location (kbp)', fontsize=size)
#plt.ylabel('coverage', fontsize=size)
plt.xlim([0,85800])
#plt.title('mtDNA coverage distribution for petite colony', fontsize=size)
#plt.grid(True)

txtsize = 20

ax.tick_params(axis='both', which='major', labelsize=size)
ax.tick_params(axis='both', which='minor', labelsize=size)
ax.tick_params(axis="y",direction="in", pad=-100)

plt.plot((82392, 82600), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(82000, text_pos, "ORI5", size = txtsize, ha='center')

plt.plot((79213, 80022), (line_pos,line_pos), 'b-', linewidth = 4.0)
plt.text(79000, text_pos_lower, "COX3", size = txtsize)

#plt.plot((74495, 75984), (line_pos,line_pos), 'r-', linewidth = 4.0)
#plt.text(74495, text_pos, "Q0255", size = txtsize)

plt.plot((73758, 74513), (line_pos,line_pos), 'b-', linewidth = 4.0)
plt.text(71000, text_pos_lower, "COX2", size = txtsize)

plt.plot((58009, 62447), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(58009, text_pos, "21S_RNA", size = txtsize)

plt.plot((56567, 56832), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(54567, text_pos_lower, "ORI4", size = txtsize)

plt.plot((54567, 54840), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(51567, text_pos, "ORI3", size = txtsize)

plt.plot((48901, 50097), (line_pos,line_pos), 'b-', linewidth = 4.0)
plt.text(47901, text_pos_lower, "VAR1", size = txtsize)

plt.plot((45227, 47927), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(45000, text_pos, "ORI6", size = txtsize)

plt.plot((36540, 43647), (line_pos,line_pos), 'b-', linewidth = 4.0)
plt.text(39540, text_pos_lower, "COB", size = txtsize)

plt.plot((32231, 32501), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(32000, text_pos, "ORI2", size = txtsize)

plt.plot((30220, 30594), (line_pos,line_pos), 'r-', linewidth = 4.0)
plt.text(29500, text_pos_lower, "ORI7", size = txtsize)

plt.plot((28487, 29266), (line_pos,line_pos), 'b-', linewidth = 4.0)
plt.text(26000, text_pos, "ATP6", size = txtsize)

#plt.plot((13818, 26701), (line_pos,line_pos), 'b-', linewidth = 4.0)
#plt.text(17000, text_pos_lower, "COX1", size = txtsize)

#plt.plot((12510, 12780), (line_pos,line_pos), 'r-', linewidth = 4.0)
#plt.text(12000, text_pos, "ORI8", size = txtsize)

#plt.plot((6546, 8194), (line_pos,line_pos), 'r-', linewidth = 4.0)
#plt.text(5700, text_pos_lower, "15S_RRNA", size = txtsize)

#plt.plot((4012, 4312), (line_pos,line_pos), 'r-', linewidth = 4.0)
#plt.text(2600, text_pos, "ORI1", size = txtsize)

fig = plt.gcf()
fig.set_size_inches(14, 9)
#plt.show()
fig.savefig(str(sys.argv[1][:-4]) + ".png", dpi=300)
print "finished: " + str(sys.argv[1][:-4]) + " in: " 
print("--- %s seconds ---" % (time.time() - start_time))
	


	
	
	
