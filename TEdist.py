#!/usr/bin/env python3

import collections
import statistics


class InternalDeletion:
    
    def __init__(self, start, end, count=0.0, freq=0.0):
        self.start = start
        self.end = end
        self.count = count
        self.freq = freq
        self.pos = (start, end)
        
    
class Truncfile:
    
    def __init__(self, filename, targetgen, targetrep):
        self.filename = filename
        self.targetgen = targetgen
        self.targetrep = targetrep
        self.allcount = collections.defaultdict(lambda:0)
        self.idcollection = collections.defaultdict(lambda:[])

    def collect_deletions_trunc(self):
        for line in open(self.filename):
            a = line.rstrip("\n").split("\t")
            rep,gen,deme,idd,start,end,count = int(a[0]),int(a[1]),int(a[2]),a[3],int(a[4]),int(a[5]),int(a[6])
            if (gen != self.targetgen or rep != self.targetrep):
                continue
            if idd == "all" and count == 0:  # disregard demes without TEs
                continue
    
            skey = deme                      # for now the key is just the deme number
            if idd == "all":
                self.allcount[skey] = count
            elif idd == "id":
                self.idcollection[skey].append(InternalDeletion(start, end, count))
            
            
class Deviatefile:
    
    def __init__(self, filename):
        self.filename = filename
        self.allcount = collections.defaultdict(lambda:0)
        self.idcollection = collections.defaultdict(lambda:[])
                
    def collect_deletions_deviate(self):
        for line in open(self.filename):
            if line.startswith('#'):
                ls = line.rstrip('\n').split(' ')
                if ls[1] == 'insertions/haploid:':
                    c = float(ls[2])
                continue
            
            a = line.rstrip('\n').split(' ')
            skey = a[1].split('/')[-1].split('.')[0] # SRR, careful!! confirmed for: NA and droseu gdl
            intdel = a[13]
            if intdel == 'NA':
                continue
            
            self.allcount[skey] = c 
            
            if ',' in intdel:
                comma_separated_deletions = intdel.split(',')
                for csd in comma_separated_deletions:
                    i = csd.split(':')
                    self.idcollection[skey].append(InternalDeletion(start=int(i[0]), end=int(i[1]), count=float(i[2])))
            else:
                i = intdel.split(':')
                self.idcollection[skey].append(InternalDeletion(start=int(i[0]), end=int(i[1]), count=float(i[2])))
                
                           
class Deme:
    
    def __init__(self, demekey, deletions, allcount):
        self.demekey = demekey
        self.allcount = allcount
        self.deletions = deletions
        
    def create_deldict(self):
        deldict = {}
        flfreq = 1.0
        
        for d in self.deletions:  # update the frequency of all internal deletions
            freq = float(d.count) / float(self.allcount)
            d.freq = freq
            # filter below a threshold frequency
            if freq > 0.005:
                flfreq -= freq
                deldict[d.pos] = freq
        
        self.deldict = deldict   # deletion dictionary
        self.flfreq = flfreq     # frequency of full length - also a marker
            
    def getFreq(self, key):
        if key in self.deldict:
            return self.deldict[key]
        else:
            return(0.0)
        
        
def heterozygosity(deletions):
    # TE = locus, deletions = alleles
    # het = 1 - sum(d_i * d_i)
    delsum = 0.0
    het = 1.0
    for d in deletions:
        if d == 0.0:
            continue
        delsum += d
        het -= d ** 2
    if delsum > 1.0:
        raise Exception("sum of deletion frequencies not allowed to exceed 1.0")
    flfreq = 1.0 - delsum
    het -= flfreq ** 2  # also include the full length as a marker
    return(het)
    

def fst(df1, df2):
    # get average of the two frequency lists
    dft = [(f1 + f2) / 2.0 for f1, f2 in zip(df1, df2)]
    h1 = heterozygosity(df1)
    h2 = heterozygosity(df2)
    ht = heterozygosity(dft)
    
    if ht == 0:
        return "na"
    hs = (h1 + h2) / 2.0
    fst = (ht - hs) / ht
    if fst < 0.0 or fst > 1.0:
        raise Exception("Invalid fst value")
    return(fst)

def unique_intdels(demes):
    # returns a unique set of intdels from all demes
    deldicts = [d.deldict for d in demes]
    idds = [list(i.keys()) for i in deldicts]
    flat_idds = [i for sublist in idds for i in sublist]
    uniq_ids = set(flat_idds)
    return(uniq_ids)
    
def construct_freqtable(demes, uniq_ids):
    freqtable = collections.defaultdict(lambda:[])
    for d in demes:
        deme_dels = list(d.deldict.keys())
        for i in uniq_ids:
            if i in deme_dels:
                freqtable[i].append(d.deldict[i])
            else:
                freqtable[i].append(0.0)
    return(freqtable)
    
def standardize(freqs):
    mean = statistics.mean(freqs)
    sd = statistics.stdev(freqs)
    st_freqs = [(x - mean) / sd for x in freqs]
    return(st_freqs)
    
def euc(df1, df2):
    euc = sum([(i - j) ** 2 for i,j in zip(df1, df2)])
    return(euc ** 0.5)

def merge(freqtable, intdel):
    candidates = list()
    for tup, freqs in freqtable.items():
        if intdel[0][0] in range(tup[0] - 3, tup[0] + 3):
            if intdel[0][1] in range(tup[1] - 3, tup[1] + 3):
                candidates.append(tup)
                
    n_cand = len(candidates)
    if n_cand > 1:
        new_start = int(round(sum([start for (start, end) in candidates]) / n_cand))
        new_end = int(round(sum([end for (start, end) in candidates]) / n_cand))
        
        cand_dict = {tup:freqs for tup, freqs in freqtable.items() if tup in candidates}
        new_freq_list = [round(sum(x),10) for x in zip(*list(cand_dict.values()))]
        return((new_start, new_end), new_freq_list, candidates)
    else:
        return(intdel[0], intdel[1], candidates)
    
def non_null(row):
    ab = len([i for i in row if i > 0.0])
    return(ab)


def calcD(df1, df2):
    
    dft = [(f1 + f2) / 2.0 for f1, f2 in zip(df1, df2)]
    h1 = heterozygosity(df1)
    h2 = heterozygosity(df2)
    ht = heterozygosity(dft)
    
    hs = (h1 + h2) / 2.0
    
    D = ( (ht - hs) / (1 - hs) ) * (2 / (2 - 1))
    if D < 0.0 or D > 1.0:
        raise Exception("Invalid D")
    return(D)
    
    
    
    