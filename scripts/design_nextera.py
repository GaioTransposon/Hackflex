#!/usr/bin/env python3
import sys
import primer3
import math

i7_file = open('len8_avoid_ill_i7_r2.txt')
i5_file = open('len8_avoid_ill_i5_r2.txt')

f7_seq = 'CAAGCAGAAGACGGCATACGAGAT'
p7_seq = 'GTCTCGTGGGCTCGG'

f5_seq = 'AATGATACGGCGACCACCGAGATCTACAC'
p5_seq = 'TCGTCGGCAGCGTC'

hairpin_tm_limit = 51
homodimer_tm_limit = 45
heterodimer_tm_limit = 45
het_i7 = {}
het_i5 = {}

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def remove_homodimers_and_hairpins(barcode_list, fseq, pseq):
    """
    Removes homodimers and hairpins from a barcode list.
    """
    good_barcodes = []
    for bc in barcode_list:
        bc_oligo = fseq+bc+pseq
        bc_hp = primer3.bindings.calcHairpin(bc_oligo,dv_conc=1,dntp_conc=0.2)
        bc_homo = primer3.bindings.calcHomodimer(bc_oligo)
        if bc_hp.tm > hairpin_tm_limit or bc_homo.tm > homodimer_tm_limit:
#            print("barcode "+bc+" has a hairpin or homodimer")
#            print(bc_hp)
#            print(bc_homo)
            continue
        else:
            good_barcodes.append(bc)
    return good_barcodes

def calc_heterodimer_network(i7_list, i5_list, quartile_tcd_tm):
    for i7 in i7_list:
        i7_oligo = f7_seq+i7+p7_seq
        het_i7[i7]={}
        for i5 in i5_list:
            if not i5 in het_i5: het_i5[i5]={}
            i5_oligo = f5_seq+i5+p5_seq
            het = primer3.bindings.calcHeterodimer(i5_oligo,i7_oligo)
            if het.tm > heterodimer_tm_limit:
                het_i7[i7][i5]=het.tm
                het_i5[i5][i7]=het.tm
#                print('over the limit')
            h2 = primer3.bindings.calcHeterodimer(i5,i7)
            if h2.tm > quartile_tcd_tm:
                het_i7[i7][i5]=h2.tm
                het_i5[i5][i7]=h2.tm
#                print('h2 over the limit')
            h3 = primer3.bindings.calcHeterodimer(i5,reverse_complement(i7))
            if h3.tm > quartile_tcd_tm:
                het_i7[i7][i5]=h3.tm
                het_i5[i5][i7]=h3.tm
#                print('h3 over the limit')

def remove_heterodimers(i7_list,i5_list,quartile_tcd_tm):
    calc_heterodimer_network(i7_list,i5_list,quartile_tcd_tm)
    new_i7=[]
    new_i5=[]
    for i7 in het_i7:
        if len(het_i7[i7])==0:
            new_i7.append(i7)
    for i5 in het_i5:
        if len(het_i5[i5])==0:
            new_i5.append(i5)
    new_lists={}
    new_lists['i5']=new_i5
    new_lists['i7']=new_i7
    return new_lists

i7_bcs = []
for line in i7_file:
    i7_bcs.append(line.rstrip())

i5_bcs = []
for line in i5_file:
    i5_bcs.append(line.rstrip())

# first screen for hairpins and homodimers
i7_bcs = remove_homodimers_and_hairpins(i7_bcs,f7_seq,p7_seq)
hairpin_tm_limit = 58
i5_bcs = remove_homodimers_and_hairpins(i5_bcs,f5_seq,p5_seq)
print(str(len(i7_bcs))+' i7 barcodes remaining after hairpin and homodimer screen')
print(str(len(i5_bcs))+' i5 barcodes remaining after hairpin and homodimer screen')

# generate a background distribution of tandem complement barcode Tm
self_tm_list=[]
for i5 in i5_bcs:
    het = primer3.bindings.calcHeterodimer(i5,reverse_complement(i5))
    self_tm_list.append(het.tm)
self_tm_list.sort()
quartile_tcd_tm = self_tm_list[int(len(self_tm_list)/4)]
print('quartile Tm: '+str(quartile_tcd_tm))

# then screen remaining oligos for heterodimers
new_lists = remove_heterodimers(i7_bcs,i5_bcs,quartile_tcd_tm)
i7_bcs = new_lists['i7']
i5_bcs = new_lists['i5']
print(str(len(i7_bcs))+' i7 barcodes remaining after hetero screen')
print(str(len(i5_bcs))+' i5 barcodes remaining after hetero screen')

wells = "ABCDEFGH"
rows = 8
cols = 12
i = 0
i7_plateplan = open('i7_plateplan.tsv','w')
for i7 in i7_bcs:
    if i >= 96: break
    well_id = wells[math.floor(i/cols)]+str(1+(i%cols))
    oligo = f7_seq+i7+p7_seq
    i7_plateplan.write(well_id+'\tuts_i7_'+reverse_complement(i7)+'\t'+oligo+'\n')
    i+=1
i7_plateplan.close()
i=0
i5_plateplan = open('i5_plateplan.tsv','w')
for i5 in i5_bcs:
    if i >= 96: break
    well_id = wells[math.floor(i/cols)]+str(1+(i%cols))
    oligo = f5_seq+i5+p5_seq
    i5_plateplan.write(well_id+'\tuts_i5_'+i5+'\t'+oligo+'\n')
    i+=1
i5_plateplan.close()
