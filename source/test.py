import numpy as np


genome_length=30
paired_read_length=3
insert_size=5
coverage=2

sequencing_fragments_num = int(round((coverage*genome_length)/(2*paired_read_length),0))

sequencing_fragments_positions = np.random.uniform(0,genome_length-insert_size,sequencing_fragments_num)



s=[round(pos) for pos in sequencing_fragments_positions]
s.sort()

print(s)


import matplotlib.pyplot as plt
count, bins, ignored = plt.hist(s, 20, density=True)
plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
plt.show()