# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:50:13 2017

@author: Owner
"""
#Functions
###############################################################################        
import random
import math
import networkx as nx

def WriteResults(list):
    """writes out elements of a list to text"""
    with open("Output.txt", "w") as text_file:
        for j in range(len(list)):
            text_file.write(str(list[j]) + " ")            
    text_file.close()
    
def readFasta(fileName):
    f = open(fileName)
    seqs = {}
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            words = line.split()
            name = words[0][1:]
            seqs[name] = ''
        else:
            seqs[name] =  seqs[name] + line
    f.close()
    return seqs
    
def readSequenceArray(fileName):
    f = open(fileName)
    seqs = []
    for line in f:
        line = line.rstrip()
        seqs.append(line)
    f.close()
    return seqs
    
def readGapSeparatedSequenceArray(fileName):
    f = open(fileName)
    seqs = []
    for line in f:
        line = line.rstrip()
        lines = line.split(' ')
        for l in lines:
            seqs.append(l)
    f.close()
    return seqs
        
def getLongestSequenceInFasta(dic):
    maxLength = 0
    longest = 0
    for k, v in dic.items():
        longest = len(v)
        if longest > maxLength:
            maxLength = longest
    return maxLength
            
def getsequenceLengthsInFasta(dic):
    sequenceL = []
    for k, v in dic.items():
        sequenceL.append(len(v))   
    return sequenceL        

def getLengthsDicInFasta(dic):
    sequenceL = {}
    for k, v in dic.items():
        l = len(v)
        sequenceL.update({k:l})   
    return sequenceL      


def reverseComplement(s):
    """computes the reverse complement of a DNA string"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def create_dna(n, alphabet='acgt'):
    """creates a random string of DNA"""
    return ''.join([random.choice(alphabet) for i in range(n)])
    
def naive(p, t):
    #Given a DNA sequence (t) search for the number of occurences of an exact 
    #match (p)
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences    
            
def kmerCount(k,t):
    #Counts the number if occurences if a k-mer in a dna sequence (naive)
    d = {}
    occurrences = []
    count = 0
    for j in range(len(t) - k + 1):
        p = t[j:j+k]
        for i in range(len(t) - len(p) + 1):        
            if(p == t[i: i + len(p)]):
                count+= 1
                occurrences.append(i) 
        d.update({p:occurrences})
        occurrences = []
    return d
    
def patternCount(t,p):
    #Counts the number if occurences if a k-mer in a dna sequence (naive)
    count = 0
    for i in range(len(t) - len(p) + 1):        
            if(p == t[i: i + len(p)]):
                count+= 1     
    return count

def approximatePatternCount(text,pattern,d):
    '''count instances of pattern with at most d mismatches'''
    count = 0
    l = []
    for i  in range(len(text) - len(pattern) + 1):
        p = text[i: i + len(pattern)]
        if(hammingDistance(pattern,p) <= d):
            count = count + 1
            l.append(i)
    return count
    
def approximatePatternMatch(text,pattern,d):
    '''count instances of pattern with at most d mismatches'''
    l = []
    for i  in range(len(text) - len(pattern) + 1):
        p = text[i: i + len(pattern)]
        if(hammingDistance(pattern,p) <= d):
           
            l.append(i)
    return l
    
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    
def getSequencingReads(filename):
    reads = []
    with open(filename, 'r') as f:
        for line in f:
            reads.append(line.rstrip())
        return reads
   
#Hamming Distance - Problem 1
def hammingDistance(s1, s2):
    """Compute the Hamming distance between two strings"""
    #test that they have values
    if len(s1) == 0: return len(s2)
    if len(s2) == 0: return len(s1)        
    #Convert the int lists to strings
    str1 = ''.join(str(e) for e in s1)
    str2 = ''.join(str(e) for e in s2)        
    #Counter set at zero
    hamDist = 0
    for i in range(0, len(str1)):
        #If the values at the specified index aren't equal
        if str1[i] != str2[i]:
            #increment
            hamDist += 1           
    #Return the total count.
    return hamDist
    
def distanceBetweenPatternAndStrings(pattern,dna):
    """Computes the the sum total distance between kmer and a list of DNA strings"""
    k = len(pattern)
    distance = 0
    totalDistance = 0
    for text in dna:
        hd = float('inf')
        # generate all kmer pattern' in text
        kmers = kmersFromDNA(text,k)
        for j in kmers:
            distance = hammingDistance(pattern,j)
            if hd > distance:
                hd = distance
        totalDistance += hd             
    return totalDistance
        
        
    
        
def medianString(dna,k):
    distance = float('inf') 
    for i in range(int(4**k)): 
        pattern = NumberToPattern(i,k)
        if distance > distanceBetweenPatternAndStrings(pattern,dna):
            distance = distanceBetweenPatternAndStrings(pattern,dna)
            median = pattern      
    return median

def LastSymbol( pattern ):
    return pattern[ -1: ]

def skew(sequence):
    """Find a position in a genome where the skew diagram attains a minimum."""
    c = 0
    g = 0
    min_skew = 0
    skew_list = []
    index = 0
    for i in sequence:
        index += 1
        if i == 'C':
            c += 1
        if i == 'G':
            g += 1
        skew = g-c
        if skew < min_skew:
            skew_list = [index]
            min_skew = skew
        if skew == min_skew and index not in skew_list:
            skew_list.append(index)    
    print(skew_list)

def SymbolToNumber( symbol ):
    retVal = 10000
    if symbol == 'A':
        retVal = 0       
    elif symbol == 'C':
        retVal = 1
    elif symbol == 'G':
        retVal = 2
    elif symbol == 'T':
        retVal = 3
    return retVal

def NumberToSymbol( index ):
    symbol = 'Z'
    if index == 0:
        symbol = 'A'
    elif index == 1:
        symbol = 'C'
    elif index == 2:
        symbol = 'G'
    elif index == 3:
        symbol = 'T'
    return symbol

def PatternToNumber( pattern ):
    """Convert a DNA string to a number"""
    if pattern == "":
        return 0
    if len( pattern ) > 0:
        subStrEndIndex = len( pattern ) - 1
    else:
        subStrEndIndex = 0
    prunedPattern = pattern[ 0: subStrEndIndex ]
    lastSymbol = LastSymbol( pattern )
    #The '4 *' allows the resulting numbers to be unique according to their symbol's positions.
    return 4 * PatternToNumber( prunedPattern ) + SymbolToNumber( lastSymbol )

def NumberToPattern( index, k ):
    """Convert an integer to its corresponding DNA string."""
    if k == 1:
        return NumberToSymbol( index )
    prefixIndex = index // 4
    remainder = index % 4
    prefixPattern = NumberToPattern( prefixIndex, k - 1 )
    symbol = NumberToSymbol( remainder )
    return prefixPattern + symbol 

def computingFrequencies(text,k):
    """Compute frequency array of kmers"""
    frequencyArray = []
    for i in range(int(4**k)):  
        frequencyArray.append(0)   
    for i in range(len(text)- k + 1):
        pattern = text[i:i+k]
        j = PatternToNumber(pattern)
        frequencyArray[j] = frequencyArray[j] + 1
    return frequencyArray
    
def frequentWords(text,k):
    """Find the most frequent k-mers in a string"""
    frequentPatterns = []
    count  = count = [0] * (len(text) - k + 1)
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        count[i] = patternCount(text,pattern)
    maxCount = max(count)
    for i in range(len(text) - k + 1):
        if(count[i] == maxCount):
            frequentPatterns.append(text[i:i+k])
    fp = set(frequentPatterns)
    result = []
    for i in fp:
        result.append(i)
    return result

def fasterFrequentWords(text,k):
    """Find the most frequent k-mers in a string"""
    frequentPatterns = []
    frequencyArray = computingFrequencies(text,k)
    maxCount = max(frequencyArray)
    for i in range(int(4**k)): 
        if frequencyArray[i] == maxCount:
            pattern = NumberToPattern(i,k)
            frequentPatterns.append(pattern)
    return frequentPatterns   

def clumpFinding(genome,k,t,L):
    """Find patterns forming clumps in a string."""
    frequentPatterns = []       
    clump = [0] * int(4**k)
    for i in range(len(genome) - L):
        text = genome[i:i + L]
        frequencyArray = computingFrequencies(text,k)
        for index in range(int(4**k)): 
            if(frequencyArray[index] >= t):
                clump[index] = 1
    for i in range(int(4**k)):
        if clump[i] == 1:
            pattern = NumberToPattern(i,k)
            frequentPatterns.append(pattern)
    return frequentPatterns
    
def neighbors(pattern,d):
    """The d-neighborhood Neighbors(Pattern, d) is the set of all k-mers whose Hamming distance from Pattern does not exceed d."""
    x = ['A','C','G','T']
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ['A','C','G','T']
    neighborhood = []
    suffixNeighbors = neighbors(suffix(pattern),d)
    for text in suffixNeighbors:
        if hammingDistance(suffix(pattern),text) < d:
            for nucleotide in x:
                p = nucleotide + text
                neighborhood.append(p)
        else:
            p = firstSymbol(pattern) + text
            neighborhood.append(p)
    return neighborhood
    
def suffix(pattern):
    return pattern[1:len(pattern)]
    
def prefix(pattern):
    return pattern[0:len(pattern) -1]

def firstSymbol(pattern):
    return pattern[0]
  
def immediateNeighbors(pattern):
    neighborhood = []
    for i in range(pattern):
        symbol = pattern[i]
        for j in range(pattern):
            if j != symbol:
                pattern[i] = j
                neighbor = pattern
                neighborhood.append(neighbor)
    return neighborhood

def frequentWordsWithMismatches(text,k,d):
    """Find the most frequent k-mers with mismatches in a string."""
    frequentPatterns = []  
    frequencyArray = []     
    close = []
    for i in range(int(4**k)):
        frequencyArray.append(0)
        close.append(0)
    for i in range(len(text) - k):
        neighborhood = neighbors(text[i:i+k],d)
        for pattern in neighborhood:
            index = PatternToNumber(pattern)
            close[index] = 1
    for i in range(int(4**k)):
        if(close[i] == 1):
            pattern = NumberToPattern(i,k)
            frequencyArray[i] = approximatePatternCount(text,pattern,d) 
    maxCount = max(frequencyArray)
    for i in range(int(4**k)):
        if frequencyArray[i] == maxCount:
            pattern = NumberToPattern(i,k)
            frequentPatterns.append(pattern)
    return frequentPatterns        
    
def motifEnumeration(dna,k,d):
    """Gets kmers from all DNA strings, searches them for kmers with at most d mistmatches then looks for d-neighbors that appear on all dna strings"""
    patterns = []
    primePatterns = []
    counter = 0
    kmers = kmersFromDNAList(dna,k) 
    # we get all neighbors of all the kmers of all dna strings
    for kmer in kmers:
        n = neighbors(kmer,d)
        for i in n:
            primePatterns.append(i)
    #now we search for all prime patterns in all string with at most d mismatches
    for j in primePatterns:
        for string in dna:
            l = approximatePatternMatch(string,j,d)
            if len(l) > 0:
                 counter += 1
        if counter >= len(dna):
            patterns.append(j) 
            counter = 0
        counter = 0   
    return set(patterns)
                                                                        
def generateKmers(k):
    """Generate all possible kmers of length k"""
    kmers = []
    for i in range(int(4**k)):
        pattern = NumberToPattern(i,k)
        kmers.append(pattern)
    return kmers

def kmersFromDNA(d,k):
    """Retunrs a list of k-mers from a string of DNA"""
    kmers = []
    for i in range(len(d) - k + 1):             
        kmers.append(d[i:i+k])        
    return kmers
    
def kmersFromDNAList(dna,k):
    kmers = []
    for i in dna:
        l = kmersFromDNA(i,k)
        for j in l:
            kmers.append(j)
    return kmers
    
def searchFrequentWords(d,l):
    max = 0
    key = ""
    value = 0
    r = {}
    for k,v in d.items():
        w = kmerCount(l,v)
        #get max
        for k,v in w.items():
            if len(v) > max:
                key = k
                value = v
                max = len(v)
        r.update({key:value})                
    return r

def genomePath(d):
    seq = d[0]
    for i in d:       
        seq += i[-1]
    return seq
       
def graphExample():
    G=nx.Graph()
    G.add_node("spam")
    G.add_edge(1,2)
    G.add_edge(2,3)
    return G
    
def kmerCompositionLexigographic(s,k):
    kmers = kmersFromDNA(s,k)
    kmers.sort()
    return kmers
    
def kmerGraph(s,k):
    d = kmersFromDNA(s,k)
    G=nx.DiGraph()
    for i in d:
        G.add_node(i)
        for j in d:
            if suffix(i) == prefix(j):
                G.add_edge(i,j)
    return G
        
###############################################################################
