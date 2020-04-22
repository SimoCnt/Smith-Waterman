"""
Es.1:
    Seq1: TFDERILGVQTYWAECLA
    Seq2: QTFWECIKGDNATY

Es.2:
    Seq1: LTGARDWEDIPLWTDWDIEQESDFKTRAFGTANCHK
    Seq2: TGIPLWTDWDLEQESDNSCNTDHYTREWGTMNAHKA
"""


import numpy as np
import Bio.SubsMat.MatrixInfo as bio
import matplotlib.pyplot as plt
import seaborn as sns
import sys


class my_dictionary(dict):

    def __init__(self):
        self = dict()

    def add(self, key, value):
        self[key] = value



################################
# Dizionario dei valori Blosum #
################################
# Contiene i valori di tutte le possibili coppie di amminoacidi per costruire la
# matrice di sostituzione
X = bio.blosum62
#X



######################################################
# Creazione matrice Blosum (matrice di sostituzione) #
######################################################
print("\nSequenza 1: ")
seq1 = input()
print("\nSequenza 2: ")
seq2 = input()
print("\nGap opening penalty: ")
delta = float(input())

bl = np.zeros((len(seq1), len(seq2)))

for i in range(bl.shape[0]):
    for j in range(bl.shape[1]):
        if(seq1[i], seq2[j]) in X:
            bl[i][j] = X[(seq1[i], seq2[j])]
        else:
            bl[i][j] = X[(seq2[j], seq1[i])]

bl = np.transpose(bl)
#bl


blNew = np.zeros((bl.shape[0] + 1, bl.shape[1] + 1))

for i in range(blNew.shape[0]):
    blNew[i][0] = 0

for j in range(blNew.shape[1]):
    blNew[0][j] = 0

for i in range(bl.shape[0]):
    for j in range(bl.shape[1]):
        blNew[i+1][j+1] = bl[i][j]

#blNew



##################################
# Creazione matrice dei punteggi #
##################################
pos = my_dictionary() # Creo un dizionario che mi conserva la provenienza di ciascun punteggio di ogni cella
mScore = np.zeros((blNew.shape[0], blNew.shape[1]))


# Celle prima riga
for j in range(mScore.shape[1]):
    mScore[0][j] = blNew[0][j]
    pos.add((0, j), ('end', mScore[0][j]))


# Celle prima colonna
for i in range(mScore.shape[0]):
    mScore[i][0] = blNew[i][0]
    pos.add((i, 0), ('end', mScore[i][0]))


# Celle centrali
for i in range(1, mScore.shape[0]):
    for j in range(1, mScore.shape[1]):
        diag = mScore[i-1][j-1] + blNew[i][j]
        sx = mScore[i][j-1] - delta
        up = mScore[i-1][j] - delta
        end = 0

        mv = max(diag, sx, up, end)
        mScore[i][j] = mv

        if(mv == diag):
            pos.add((i,j), ('diag', mv))
        elif(mv == sx):
            pos.add((i,j), ('sx', mv))
        elif(mv == up):
            pos.add((i,j), ('up', mv))
        else:
            pos.add((i,j), ('end', mv))

#pos


# Creo una sottomatrice subMScore di mScore (mScore meno prima riga meno prima colonna)
subMScore = mScore[1:mScore.shape[0], 1:mScore.shape[1]]


# Inizializzo i due indici alla posizione del punteggio massimo della matrice mScore
max = 0
for i in range(subMScore.shape[0]):
    for j in range(subMScore.shape[1]):
        if(subMScore[i][j] >= max):
            max = subMScore[i][j]
            i_max = i
            j_max = j


# Punteggio dell'allineamento
punteggio = subMScore[i_max][j_max]


"""
plt.figure(figsize=(9,9))
sns.heatmap(mScore, annot=True, cmap="Blues_r", linewidths=.5, square=True, xticklabels=seq1, yticklabels=seq2)
"""



#############################
# Ricerca percorso migliore #
#############################
percorso = []

i = i_max
j = j_max

ax = []
ay = []


if(pos[(i, j)][0] == 'diag'):
    ax.append(i+1)
    ay.append(j+1)
elif(pos[(i, j)][0] == 'sx'):
    ax.append(i)
    ay.append(j+1)
else:
    ax.append(i+1)
    ay.append(j)

while True:
    if(pos[(i, j)][0] == 'diag'):
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        i = i - 1
        j = j - 1
    elif(pos[(i, j)][0] == 'sx'):
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        j = j - 1
    elif(pos[(i, j)][0] == 'up'):
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        i = i - 1
    else: # Termina quando pos[(i, j)][0] == 'end'
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        break

percorso.reverse()
#percorso



################################
# Allineamento tra le stringhe #
################################
newSeq1 = []
newSeq2 = []
confr = []

percorso.reverse()
percorso

for i in range(len(percorso)):
    if(percorso[i][1] == 'diag'):
        newSeq1.append(percorso[i][0][0])
        newSeq2.append(percorso[i][0][1])
    elif(percorso[i][1] == 'up'):
        newSeq1.append('-')
        newSeq2.append(percorso[i][0][1])
    elif(percorso[i][1] == 'sx'):
        newSeq1.append(percorso[i][0][0])
        newSeq2.append('-')
    else:
        newSeq1.append(percorso[i][0][0])
        newSeq2.append(percorso[i][0][1])
        break

newSeq1.reverse()
newSeq2.reverse()


if(len(newSeq1) < len(newSeq2)):
    count = len(newSeq1)
else:
    count = len(newSeq2)

for i in range(count):
    if(newSeq1[i] == newSeq2[i] and newSeq1[i] != '-' and newSeq2 != '-'):
        confr.append('|')
    else:
        confr.append(' ')



#####################################
# Stampa dell'allineamento migliore #
#####################################
print('\n')
print(*newSeq1, sep=" ")
print(*confr, sep=" ")
print(*newSeq2, sep=" ")
print('\n')



#####################
# Plot del Percorso #
#####################
#punteggio = round(punteggio)

plt.figure(figsize=(12,12))
plt.title("Percorso migliore (punteggio: {})".format(punteggio))
sns.heatmap(subMScore, annot=True, cmap="Blues_r", linewidths=.5, square=True, xticklabels=seq1, yticklabels=seq2)
plt.plot(ay, ax, 'r-')
plt.show()
