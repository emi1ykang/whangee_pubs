# Import all needed modules
import cPickle as pic
import numpy as np
import matplotlib.pyplot as plt

# Setup pickles to be imported
allele = pic.load(open("allele_dic.pkl", "rb"))
number = pic.load(open("aminotonumber.pkl", "rb"))
translate = pic.load(open("translate.pkl", "rb"))

# Array of appropiate size, filled with zeros
heat1 = np.zeros((21, 78))

# Pull & Sort data from pickle files
# for loop iterate over each mutant that was sequenced
for key in allele:
    
    # residue that was mutated and what it was mutated to in DNA space
    # deconstruct from pickle into a tuple
    (pos,codon) = (allele[key][0]).split('_')
    
    # 'transcribe' mutation into RNA space
    codon = codon.replace('T', 'U')
    
    # 'translate' mutation into protein space
    aa = translate[codon]
    
    # link each amino acid and stop codon to a number 1-21
    # this makes plotting in an array easier
    assign = number[aa]
    
    # change residue position from a string to an integer
    pos = int(pos)
    
    # positions were recorded as 2-77, and it is easier to begin plotting at 0
    # note: ubiquitin has 76 residues...label accordingly
    pos = pos - 2
    
    # fill array according to frequency at each position in space
    heat1[assign][pos]+=1

# Assigns color according to values in heat1 array
heatmap1 = plt.pcolor(heat1)

# Sorting out labels for y-axis in correct order and placing in new list 'yax'
yax1 = sorted(number.items(), key=lambda (k, v): v)
yax = []
for item in yax1:
    yax.append (item[0])

# Setting up y-axis ticks...want AA labels, centered
y = np.arange(0.5, 21.5, 1)
plt.yticks(y, yax, rotation='horizontal')

# Setting up x-axis ticks...want residue number, centered
x = np.arange(0.5, 76.5, 5)
xlab = np.arange(1, 77, 5)
plt.xticks(x, xlab, rotation='horizontal')

# Limit axes to appropriate values
plt.axis([0, 76, 0, 21])

plt.show()
