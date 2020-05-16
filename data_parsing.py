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