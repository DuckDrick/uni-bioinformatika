CODONS = {
    'start': ['ATG'],
    'stop': ['TAA', 'TAG', 'TGA'],
}

from math import sqrt
from Bio import SeqIO;
import re;
import numpy

def findSubstringLocation(a_str, sub):
    locations = [];
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            break;
        locations.append(start);
        start += len(sub);
    return locations;

def groupByOffset(locations):
    zero = filter(lambda l: l % 3 == 0 ,locations);
    one = filter(lambda l: l % 3 == 1 ,locations);
    two = filter(lambda l: l % 3 == 2 ,locations);
    return [list(zero), list(one), list(two)];

def findCodonPairs(sequence):
    starts = [];
    
    for codon in CODONS['start']:
        starts.append(findSubstringLocation(sequence, codon))


    stops = [];
    
    for codon in CODONS['stop']:
        stops.append(findSubstringLocation(sequence, codon))

    starts = [item for sublist in starts for item in sublist]
    stops = [item for sublist in stops for item in sublist]

    starts.sort();
    stops.sort();

    starts = groupByOffset(starts)
    stops = groupByOffset(stops)


    rezults = [];
    for offset in range(3):
        rez = {}
        startIndex = 0;
        stopIndex = 0;

        startPositions = starts[offset];
        stopPositions = stops[offset];
        while stopIndex < len(stopPositions) and startIndex < len(startPositions):
            if stopPositions[stopIndex] > startPositions[startIndex]:
                if stopPositions[stopIndex] not in rez:
                    rez[stopPositions[stopIndex]] = [startPositions[startIndex]]
                else:
                    rez[stopPositions[stopIndex]].append(startPositions[startIndex])
                startIndex += 1;
            else:
                stopIndex += 1;
        rezults.append(rez)

    findFurthestCodonPairs(rezults)
    return rezults;

def findFurthestCodonPairs(positions):
    rezults = [];
    for group in positions:
        rez = []
        for stopPositions in group:
            rez.append([min(group[stopPositions]), stopPositions]);
        rezults.append(rez);

    filterOutShorterThan(rezults, 100)
    return rezults;
    
def filterOutShorterThan(pairs, length):
    rezults = []
    for group in pairs:
        rezults.append(list(filter(lambda pair: pair[1] - pair[0] + 3 >= length, group)))

    frequencies(rezults)
    return rezults;

codonFrequencies = {};
didoconFrequencies = {};
aminoFrequencies = {};
diaminoFrequencies = {};

def frequencies(pairs):
    for group in pairs:
        for pair in group:
            for offset in range(1):
                for codon in re.findall(r".{1,3}", str(sequence[pair[0] + 3 + offset :pair[1]])):
                    if len(codon) == 3:
                        codonFrequencies[codon] += 1
                for dicodon in re.findall(r".{1,6}", str(sequence[pair[0] + 3 + offset :pair[1]])):
                    if len(dicodon) == 6:
                        didoconFrequencies[dicodon] += 1

                for amino in str(sequence[pair[0] + 3 + offset :pair[1]].translate()):
                    aminoFrequencies[amino] += 1
                for amino in re.findall(r".{1,2}", str(sequence[pair[0] + 3 + offset :pair[1]].translate())):
                    if (len(amino) == 2):
                        diaminoFrequencies[amino] += 1

sequence = '';



def genAllCodons():
    vals = ['A', 'T', 'G', 'C']
    for first in vals:
        for second in vals:
            for third in vals:
                codonFrequencies[first + second + third] = 0;
                for fourth in vals:
                    for fifth in vals:
                        for sixth in vals:
                            didoconFrequencies[first + second + third + fourth + fifth + sixth] = 0;
    aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*'];
    for first in aminos:
        aminoFrequencies[first] = 0
        for second in aminos:
            diaminoFrequencies[first+second] = 0;





failai = ['./data/bacterial1.fasta', './data/bacterial2.fasta', './data/bacterial3.fasta', './data/bacterial4.fasta',
'./data/mamalian1.fasta', './data/mamalian2.fasta', './data/mamalian3.fasta', './data/mamalian4.fasta',]
genAllCodons();

headerscodon = ['']
for freq in codonFrequencies:
    headerscodon.append(freq);

headersdicodon = ['']
for freq in didoconFrequencies:
    headersdicodon.append(freq);

headersamino = ['']
for freq in aminoFrequencies:
    headersamino.append(freq);

headersdiamino = ['']
for freq in diaminoFrequencies:
    headersdiamino.append(freq);

dataCodon = [headerscodon];
dataDicodon = [headersdicodon];
dataAmino = [headersamino];
dataDiamino = [headersdiamino];

for failas in failai:
  
    genAllCodons();
    fasta_sequences = SeqIO.parse(open(failas),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq;
        findCodonPairs(sequence);

        sequence = sequence.reverse_complement();
        findCodonPairs(sequence);

        totalcodons = 0;
        totaldicodons = 0;
        totalaminos = 0;
        totaldiaminos = 0;
        for freq in codonFrequencies:
            totalcodons += codonFrequencies[freq];
        for freq in didoconFrequencies:
            totaldicodons += didoconFrequencies[freq];
        for freq in aminoFrequencies:
            totalaminos += aminoFrequencies[freq];
        for freq in diaminoFrequencies:
            totaldiaminos += diaminoFrequencies[freq];


        for freq in codonFrequencies:
            codonFrequencies[freq] /= totalcodons;
    
        for freq in didoconFrequencies:
            didoconFrequencies[freq] /= totaldicodons;
        
        for freq in aminoFrequencies:
            aminoFrequencies[freq] /= totalaminos;

        for freq in diaminoFrequencies:
            diaminoFrequencies[freq] /= totaldiaminos;

        kazkas = []
        kazkas.append(failas.split('/')[2])
        for freq in codonFrequencies:
            kazkas.append(codonFrequencies[freq]);

        dataCodon.append(kazkas);
      
        kazkas = []
        kazkas.append(failas.split('/')[2])
        for freq in didoconFrequencies:
            kazkas.append(didoconFrequencies[freq]);

        dataDicodon.append(kazkas)

        kazkas = []
        kazkas.append(failas.split('/')[2])
        for freq in aminoFrequencies:
            kazkas.append(aminoFrequencies[freq]);

        dataAmino.append(kazkas)

        kazkas = []
        kazkas.append(failas.split('/')[2])
        for freq in diaminoFrequencies:
            kazkas.append(diaminoFrequencies[freq]);

        dataDiamino.append(kazkas)


t_matrix = numpy.transpose(dataCodon)

f = open("codonFrequencies.csv", "w")
for line in t_matrix:

    f.write(','.join(line));
    f.write('\n');
f.close()


t_matrix = numpy.transpose(dataDicodon)
f = open("dicodonFrequencies.csv", "w")
for line in t_matrix:
    f.write(','.join(line));
    f.write('\n');
f.close()

t_matrix = numpy.transpose(dataAmino)
f = open("aminoFrequencies.csv", "w")
for line in t_matrix:
    f.write(','.join(line));
    f.write('\n');
f.close()

t_matrix = numpy.transpose(dataDiamino)
f = open("diaminoFrequencies.csv", "w")
for line in t_matrix:
    f.write(','.join(line));
    f.write('\n');
f.close()

print(len(dataCodon) - 1)
for i in range(1, len(dataCodon)):
    print(dataCodon[i][0], end='');
    for j in range(1, len(dataCodon)):
        nzn = 0
        for k in range(1, len(dataCodon[i])):
            nzn += (dataCodon[i][k] - dataCodon[j][k]) ** 2;

        nzn = sqrt(nzn);
        print(f' {nzn}', end='');
    print('')

print('\n')
print(len(dataDicodon) - 1)
for i in range(1, len(dataDicodon)):
    print(dataDicodon[i][0], end='');
    for j in range(1, len(dataDicodon)):
        nzn = 0
        for k in range(1, len(dataDicodon[i])):
            nzn += (dataDicodon[i][k] - dataDicodon[j][k]) ** 2;

        nzn = sqrt(nzn);
        print(f' {nzn}', end='');
    print('')
print('\n')

print(len(dataAmino) - 1)
for i in range(1, len(dataAmino)):
    print(dataAmino[i][0], end='');
    for j in range(1, len(dataAmino)):
        nzn = 0
        for k in range(1, len(dataAmino[i])):
            nzn += (dataAmino[i][k] - dataAmino[j][k]) ** 2;

        nzn = sqrt(nzn);
        print(f' {nzn}', end='');
    print('')
print('\n')

print(len(dataDiamino) - 1)
for i in range(1, len(dataDiamino)):
    print(dataDiamino[i][0], end='');
    for j in range(1, len(dataDiamino)):
        nzn = 0
        for k in range(1, len(dataDiamino[i])):
            nzn += (dataDiamino[i][k] - dataDiamino[j][k]) ** 2;

        nzn = sqrt(nzn);
        print(f' {nzn}', end='');
    print('')
