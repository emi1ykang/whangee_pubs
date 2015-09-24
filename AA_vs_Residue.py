# Import all needed modules
import cPickle as pic
import numpy as np
import matplotlib.pyplot as plt

# Setup pickles to be imported
allele = pic.load(open("allele_dic.pkl", "rb"))
number = pic.load(open("aminotonumber.pkl", "rb"))
translate = pic.load(open("translate.pkl", "rb"))

# Array filled with zeros
heat1 = np.zeros((21, 78))

# Pull & Sort data from pickle files
for key in allele:
    (pos,codon) = (allele[key][0]).split('_')
    codon = codon.replace('T', 'U')
    
    aa = translate[codon]
    
    assign = number[aa]
    
    pos = int(pos)
    
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

# Setting up x-axis ticks...want residue #, centered
x = np.arange(0.5, 77.5, 5)
xlab = np.arange(1, 78, 5)
plt.xticks(x, xlab, rotation='horizontal')

# Limit axes to appropriate values
plt.axis([0, 78, 0, 21])

plt.show()