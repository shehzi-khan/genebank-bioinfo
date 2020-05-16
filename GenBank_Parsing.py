"http://biopython.org/DIST/docs/tutorial/Tutorial.html"
import os
import pandas
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import csv

keys_match={
        "Pid":"protein_id",
        "Name":"gene",
        "Name_alt":"locus_tag",
        "Product":"product",
        "Note":"note",
        "Seqa":"translation",
    }

orf_file = pandas.read_csv("CP001509_170113.orf",sep='\t')

for row in orf_file.iterrows():
    print(row)
    break

orfs=[]
keys=[]
record = SeqIO.read("CP001509_20180315.gbk.txt", "genbank")
for i,feature in enumerate(record.features):
    if feature.type!="source":
        matches=[doc for doc in orfs if doc["Begin"]==int(feature.location.start) and doc["End"]==int(feature.location.end)]

        if len(matches)==0:
            doc={
                "Genome_accn":record.id,
                "Id_orf":i+1,
                "Begin":int(feature.location.start),
                "End":int(feature.location.end),
                "Strand":"+" if feature.location.strand>0 else "-",
                "Size":int(feature.location.end)-int(feature.location.start),
                "Seqn":record.seq[int(feature.location.start):int(feature.location.end)]
            }

            orfs.append(doc)
        else:

            match=matches[0]
            for key in keys_match.keys():
                # print(key, feature.qualifiers.keys())
                if keys_match[key] in feature.qualifiers.keys():

                    match[key]=feature.qualifiers[keys_match[key]]

for rc in orfs:
    print(rc)
# output_file=open("output_file.orf", 'w')
# writer = csv.DictWriter(output_file, fieldnames=orfs[0].keys())
# writer.writeheader()
# for data in orfs:
#     writer.writerow(data)
exit()



all_types=[feature.type for feature in record.features]

cds=[]
genes=[]
sources=[]
other=[]
total=0
feature_types=set(all_types)
for type in feature_types:
    print(type,all_types.count(type))
    if type!="gene":
        total=total+all_types.count(type)
print(total)
exit()
for feature in record.features:
    if feature.type=="CDS":
        cds.append(feature)
    elif feature.type=="gene":
        genes.append(feature)
    elif feature.type=="source":
        sources.append(feature)
    else:
        other.append(feature)


print(len(cds), len(genes),len(sources),len(other),orf_file.size,orf_file.shape)
for feature in other:
    print(feature.type)
exit()
print(record.features[1])
table = 11
min_pro_len = 100

exit()
def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((record.id,start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    answer.sort()
    return answer

orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
print(orf_list[0])