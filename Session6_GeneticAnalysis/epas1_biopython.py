# coding: utf-8

# Biopython tutorial and sequence objects

# Written by R. Antonio Herrera for AGAR - 2018 Workshop

print('importing packages\n')
import numpy as np
import pandas as pd
import time

from Bio import motifs
from Bio.Seq import Seq
from Bio import SearchIO
from Bio import ExPASy
from Bio import SwissProt


### To access Entrez databases:
from Bio import Entrez
Entrez.email = "your@email.com"
handle = Entrez.einfo()
record =  Entrez.read(handle)
record.keys()
print('list of databases in Entrez:')
for l in record['DbList']:
    print(l)

print('sleeping to slow down calls to Entrez')
zs = 'z'
nap_time = np.arange(0,6)
for period in nap_time:
    print(zs)
    zs = str(zs) + 'z'
    time.sleep(5)

### Fetching gene information using a Gene ID
from Bio import Entrez
import time
Entrez.email ="eigtw59tyjrt403@gmail.com"

# epas1 ID
epas1_id = 'U81984.1'

# upload a list of IDs beforehand
gis=[epas1_id]

# form a request and choose the database
request = Entrez.epost("nucleotide",id=",".join(map(str,gis)))

# fetch the results
result = Entrez.read(request)
webEnv = result["WebEnv"]
queryKey = result["QueryKey"]

# request the results in xml format
handle = Entrez.efetch(db="nucleotide",retmode="xml", webenv=webEnv, query_key=queryKey)

# parse the results using biopython
results = []
for r in Entrez.parse(handle):
    # Grab the GI
    try:
        gi=int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
    except ValueError:
        gi=None
    print("\n>GI ",gi," "+r["GBSeq_primary-accession"]+" "+r["GBSeq_definition"]+"\n")
    results.append(r)
for key in r.keys():
    print(key)


print('First 60 nucleotides of EPAS1:\n', results[0]['GBSeq_sequence'][0:60])

epas1_dna = r['GBSeq_sequence']

print('full EPAS1 sequence\n', epas1_dna)

### Do not perform these operations during peak NCBI hours, BLAST can be run locally too
# from Bio.Blast import NCBIWWW
# from Bio.Blast import NCBIXML

# result_handle = NCBIWWW.qblast("blastn", "nt", str(epas1_id))

# with open(str(epas1_id) + "_blast.xml", "w") as out_handle:
#     out_handle.write(result_handle.read())

# result_handle.close()

# result_handle = open(str(epas1_id) + "_blast.xml")

# blast_records = NCBIXML.read(result_handle)
# E_VALUE_THRESH = 0.04

# for alignment in blast_records.alignments:
#     for hsp in alignment.hsps:
#         if hsp.expect < E_VALUE_THRESH:
#             print('****Alignment****')
#             print('sequence:', alignment.title)
#             print('length:', alignment.length)
#             print('e value:', hsp.expect)
#             print(hsp.query[0:75] + '...')
#             print(hsp.match[0:75] + '...')
#             print(hsp.sbjct[0:75] + '...')

# blast_qresult = SearchIO.read(str(epas1_id) + "_blast.xml", 'blast-xml')

# # print(blast_qresult)
# blast_hit = blast_qresult[3]
# print(blast_hit)

# suppose we have sequences for motifs of a given transcription factor binding site
instances = [Seq("TACAA"), Seq("TACGC"), Seq("TACAC"), Seq("TACCC"), Seq("AACCC"), Seq("AATGC"), Seq("AATGC")]

# using the motifs module create m motifs of the instances of the sequences
m = motifs.create(instances)

# print all motifs
print('motifs')
print(m)

# nucleotide counts can be accessed on whole
print('\nnucleotide counts')
print(m.counts)

# or by nucleotide, by using the base name as the dict key
print('\nonly counts of A')
print(m.counts['A'])

# or access all elements in the first two positions of the motifs
print('\nfirst two positions in motif collection')
print(m.counts[:,0:1])

# or print the consensus
print('\nconsensus sequence')
print(m.consensus)

# or the anticonsensus
print('\nanticonsensus sequence')
print(m.anticonsensus)

m.weblogo("my_motif.png",kwargs=['stack_width:medium','stacks_per_line:40','alphabet:alphabet_dna','ignore_lower_case:True','unit_name:bits','first_index:1',
'logo_start:1',
'composition:comp_auto',
'show_xaxis:False',
'show_yaxis:False',
'yaxis_scale:auto',
'yaxis_tic_interval:1.0',
'color_scheme:color_auto',])

print('\nweb logo motif saved as my_motif.png')

accessions = ["Q99814"]
for accession in accessions:
    handle = ExPASy.get_sprot_raw(accession)
    record = SwissProt.read(handle)
    handle = ExPASy.get_sprot_raw(accession)
    text= handle.read()
data = [x.strip() for x in text.split('\n')]
aa = data[-17:]
AA = ''
for a in aa:
    AA = str(AA) + str(a)
AA = AA.split(' ')
new_aa = ''.join(AA)[0:-2]

print('Epas1 amino acid sequence:\n'+str(new_aa))
