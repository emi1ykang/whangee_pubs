import cPickle as pic
import numpy as np
import matplotlib.pyplot as plt
import seaborn

alleledict = pic.load(open("allele_dic.pkl"))
translatedict = pic.load(open("translate.pkl"))
aa2numdict = pic.load(open("aminotonumber.pkl"))

print alleledict
print translatedict

codontonumdict = {}
index = 0
for codon in translatedict:
	codontonumdict[codon] = index
	index += 1

'''codondict = {}
for z in range(0,64):'''

posdict = np.zeros((78, 21))
codondict = np.zeros((78, 64))

for x in alleledict:
	(pos, codon) = alleledict[x][0].split('_')
	codon = codon.replace('T', 'U')
	aminoacid = translatedict[codon]
	aanum = aa2numdict[aminoacid]
	posdict[int(pos)][int(aanum)] += 1
	codondict[int(pos)][codontonumdict[codon]] += 1

#sumplot = codondict.sum(axis=0)
heatmaptop = plt.pcolor(posdict)
heatmaptop2 = plt.pcolor(codondict)
#heatmaptop3 = plt.hist(sumplot)

#plt.show(heatmaptop)
plt.show(heatmaptop2)
#plt.show(heatmaptop3)