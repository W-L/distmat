#!/usr/bin/env python3

import argparse
import TEdist
import pandas
import collections
    
    
parser = argparse.ArgumentParser()
parser.add_argument("--id-file", type=str, required=True, dest="idfile", help="Internal deletion file")
parser.add_argument("--filetype", type=str, required=True, dest="filetype", help="either invade or deviate")
parser.add_argument("--replicate", type=int, dest="rep", default=1)
parser.add_argument("--generation", type=int, dest="gen", default=1)
parser.add_argument("--distance-measure", type=str, dest="distm", default="D")
args = parser.parse_args()


if args.filetype == 'invade':
    sample = TEdist.Truncfile(filename=args.idfile, targetgen=int(args.gen), targetrep=int(args.rep))
    sample.collect_deletions_trunc() # fill the allcount and idcollection dicts
elif args.filetype == 'deviate':
    sample = TEdist.Deviatefile(filename=args.idfile)
    sample.collect_deletions_deviate() # fill the allcount and idcollection dicts
else:
    raise Exception('invalid filetype')

# get demes
demes = []
for key, count in sample.allcount.items():
    d = TEdist.Deme(demekey=key, allcount=count, deletions=sample.idcollection[key])
    d.create_deldict()
    demes.append(d)
demes.sort(key=lambda d:d.demekey)


# write out the TE abundance and FL-freq
if args.filetype == 'deviate':
    flfreqs = open(args.idfile + '.flfreq', 'w')
    flfreqs.write('pop' + ' ' + 'flfreq' + ' ' + 'copyn' + '\n')
    for d in demes:
        flfreqs.write(str(d.demekey) + ' ' + str(d.flfreq) + ' ' + str(d.allcount) + '\n')
flfreqs.close()


# create frequency table, key=deletion id, val=list of frequencies
unique_ids = TEdist.unique_intdels(demes)
freqtable = TEdist.construct_freqtable(demes, unique_ids)

# merge similar deletions
merged_freqtable = {}
while len(freqtable) != 0:
    intdel = next(iter(freqtable.items()))
    pos, freq, candidates = TEdist.merge(freqtable, intdel)
    for k in candidates:
        freqtable.pop(k, None)
    
    if pos not in merged_freqtable.keys():
        # filter deletions that appear in >n demes
        if TEdist.non_null(freq) > 1:
            merged_freqtable[pos] = freq
    else:
        if merged_freqtable[pos] != freq:
            print(merged_freqtable[pos])
            print(freq)
            
            

# save the frequency table
cols = [d.demekey for d in demes]
freqtable_frame = pandas.DataFrame.from_dict(merged_freqtable, orient='index', columns=cols)
freqtable_frame.index.name = 'deletion'
freqtable_frame.to_csv(args.idfile + '.freqtable', sep=' ', mode='w')


# calc DFS
counts = []
for intdel, freqs in merged_freqtable.items():
    f = TEdist.non_null(freqs)
    counts.append(f)
    
dfs = open(args.idfile + '.dfs', 'w')
dfs.write('ab' + ' ' + 'fre' + '\n')
for ab, fre in collections.Counter(counts).items():
    dfs.write(str(ab) + ' ' + str(fre) + '\n')
dfs.close()

# standardize deletion frequencies
st_freqtable = {key:TEdist.standardize(val) for (key,val) in merged_freqtable.items()}
# add std freqs to demes
for i in range(0, len(demes)):
    st_f = [val[i] for (key,val) in st_freqtable.items()]
    f = [val[i] for (key,val) in merged_freqtable.items()]
    demes[i].st_freqs = st_f
    demes[i].freqs = f
    


# compute pairwise distance matrix; initialize with zero - list of lists
dm = [[0.0, ] * len(demes) for i in range(0, len(demes))]

for i in range(0, len(demes)):
    for k in range(0, len(demes)):
        if args.distm == 'fst':
            dist_fst = TEdist.fst(demes[i].freqs, demes[k].freqs)
            dm[i][k] = dist_fst
        elif args.distm == 'euc':
            dist_euc = TEdist.euc(demes[i].st_freqs, demes[k].st_freqs)
            dm[i][k] = dist_euc
        elif args.distm == 'D':
            dist_D = TEdist.calcD(demes[i].freqs, demes[k].freqs)
            dm[i][k] = dist_D
        

# open file object
if args.filetype == 'invade':
    file_distm = open(args.idfile + '_rep' + str(args.rep) + '_gen' + str(args.gen) + '.' + str(args.distm), 'w')
elif args.filetype == 'deviate':
    file_distm = open(args.idfile + '.' + str(args.distm), 'w')

# write a header
header = [''] + [str(d.demekey) for d in demes]
file_distm.write('\t'.join(header) + '\n')

# write the distances
for i in range(0, len(demes)):    
    li = [demes[i].demekey] + dm[i]
    line = [str(k) for k in li]
    file_distm.write('\t'.join(line) + '\n')
    
file_distm.close()

    



