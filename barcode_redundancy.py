### Barcode Redundancy Chart...works!
# Import all needed modules
import cPickle as pic
import numpy as np
import matplotlib.pyplot as plt

from collections import Counter

# Setup pickles to be imported
allele = pic.load(open("allele_dic.pkl", "rb"))
number = pic.load(open("aminotonumber.pkl", "rb"))
translate = pic.load(open("translate.pkl", "rb"))


codons = []
# Pull & Sort data from pickle files
# for loop iterate over each mutant that was sequenced
for key in allele:
    
    # residue that was mutated and what it was mutated to in DNA space
    # deconstruct from pickle into a tuple
    (pos, codon) = (allele[key][0]).split('_')
    
    # 'transcribe' mutation into RNA space and store in new list
    codon = codon.replace('T', 'U')
    codons.append(codon)

# Count occurences of each codon, and sort alphabetically
# Put into a list of tuples...eg: [('AAA', 653), ('AAC, 400), ....]
c = Counter( codons )
sorted_by_codon = sorted(c.items(), key=lambda tup: tup[0])

# Set up bar chart to plot number of barcodes mapping to each codon
y = []
for item in sorted_by_codon:
    y.append(item[1])
    
x = np.arange(0, 64, 1)

plt.bar(x, y)

# Setting up x-axis ticks...want codon, centered
xlab = []
for item in sorted_by_codon:
    xlab.append(item[0])
xax = np.arange(0.5, 64.5, 1)
plt.xticks(xax, xlab, rotation='vertical')

# Setting up y-axis ticks
yax = np.arange(0, 1201, 100)
plt.yticks(yax, yax, rotation='horizontal')

# Other chart parameters
plt.axis([0, 64, 0, 1200])
plt.xlabel("Codon")
plt.ylabel("Number of Bar Codes")
plt.title("Redundancy of Bar Codes")
plt.grid(True)
plt.show()