#!/usr/bin/env python3
import argparse

import sys

def getOptions(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description="Parses command.")

    parser.add_argument("-i", "--input", required=True, help="Your input fasta file.")

    parser.add_argument("-o", "--output", required=True, help="Your output path.")

    parser.add_argument("-p", "--prefix",default='ORF', help="Your output name prefix.")

    options = parser.parse_args(args)

    return options


options = getOptions(sys.argv[1:])

import pyfastx

import sys

import re

import pandas as pd

d = """TTT F      CTT L      ATT I      GTT V
TTC F      CTC L      ATC I      GTC V
TTA L      CTA L      ATA I      GTA V
TTG L      CTG L      ATG M      GTG V
TCT S      CCT P      ACT T      GCT A
TCC S      CCC P      ACC T      GCC A
TCA S      CCA P      ACA T      GCA A
TCG S      CCG P      ACG T      GCG A
TAT Y      CAT H      AAT N      GAT D
TAC Y      CAC H      AAC N      GAC D
TAA Stop   CAA Q      AAA K      GAA E
TAG Stop   CAG Q      AAG K      GAG E
TGT C      CGT R      AGT S      GGT G
TGC C      CGC R      AGC S      GGC G
TGA Stop   CGA R      AGA R      GGA G
TGG W      CGG R      AGG R      GGG G"""

d = dict(zip(d.split()[::2],d.split()[1::2]))

#ini_codon=['ATG','CTG','GTG','TTG','ACG','AAG','AGG','ATC','ATT','ATA']
#ini_codon=['ATG','CTG','GTG','TTG']
ini_codon=['ATG']

########## 通过起始密码子与终止密码子切割
def cut(rna,ini_codon):
  orflist=[];
  for ic in ini_codon:
    if ic in rna:
      for m in re.finditer(ic,rna):
        seqlist=rna[m.start()::];
        seqlist = re.findall(r'.{3}',seqlist);
        orf='';
        if 'TAA' in seqlist or 'TAG' in seqlist or 'TGA' in seqlist:
          for i in seqlist:
            orf=orf+i;
            if i in ['TAA','TAG','TGA']:
              break;
          orflist.append(orf);
  return orflist;


########翻译ORF
def translate(dna):
  aa='';
  for i in range(0, len(dna), 3):
    codon=dna[i:i+3];
    if d[codon]=='Stop':
      break;
    aa += d[dna[i:i+3]];
  return aa;


#########生成fasta文件

def orf_find(fastafile=options.input, outpath=options.output,prefix=options.prefix):
  with open(outpath+"/"+prefix+".fa","w") as orf_out:
    for name,seq_rna in pyfastx.Fasta(fastafile,build_index=False):
      orflist=cut(rna=seq_rna,ini_codon=ini_codon);
      if len(orflist) == 0:
        print('transcript '+name+' does not have orf');
        continue;
      # Preserve discovery order so repeated runs produce the same FASTA order.
      orflist=list(dict.fromkeys(orflist));
      orflist=[ x for x in orflist if len(x)>=18 ];
      if len(orflist) == 0:
        print('transcript '+name+' does not have orf which length over 5!');
        continue;
      peplist=[];
      for i in range(len(orflist)):
        if len(re.findall('A',orflist[i]))+len(re.findall('T',orflist[i]))+len(re.findall('C',orflist[i]))+len(re.findall('G',orflist[i]))==len(orflist[i]):
          for m in re.finditer(orflist[i],seq_rna):
            orf_start=m.start();
            orf_end=m.end();
            orf_len=(orf_end - orf_start);
            seq_pep=translate(orflist[i]);
            peplist.append([orf_start,orf_end,orf_len,seq_pep]);
      if len(peplist) == 0:
        print('transcript '+name+' does not have an ORF containing only A/C/G/T bases');
        continue;
      peplist=pd.DataFrame(peplist,columns=('start','end','length','seq'));
      peplist_simple=peplist.groupby('end').apply(lambda x:x[x['length']==max(x['length'])]);
      peplist_simple['start']=peplist_simple['start'].astype('str');
      peplist_simple['end']=peplist_simple['end'].astype('str');
      peplist_simple.apply(lambda x:orf_out.write(">lcl|"+name+":"+x['start']+":"+x['end']+"\n"+x['seq']+"\n"),axis=1);
  orf_out.close();


if __name__ == '__main__':

  orf_find()
