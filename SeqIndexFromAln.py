#!/usr/bin/env python
import sys
(fasta_file, ref_seq, index_file) = sys.argv[1:]
def read_fasta (fp):
    seq_id, name, seq = None, None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))

fp = open (fasta_file, 'r')
alignment, reference  = {}, {}
for seq_id, seq in read_fasta (fp):
    name = seq_id[1:]
    #print name, seq
    if name == ref_seq:
        aln_num, ref_seq_num = 0, 0
        for i in seq:
            aln_num += 1
            if i != '-':
                ref_seq_num += 1
            reference[aln_num] = ref_seq_num
            #print name, aln_num, ref_seq_num
    else:
        alignment[name] = {}
        aln_num, seq_num = 0, 0
        for i in seq:
            aln_num += 1
            if i != '-':
                seq_num += 1
                alignment[name][aln_num] =  seq_num
                #print name, aln_num, alignment[name][aln_num]
with open(index_file, "w") as out_file:
    out_file.write ("\t".join(('pdb', 'pdb_seq_num', 'ref_seq_num')))
    out_file.write("%s" % "\n")
    for pdb in alignment:
        for aln_num in sorted(alignment[pdb], key = int):
            if aln_num in reference:
                ref_num = reference[aln_num]
                seq_num = alignment[pdb][aln_num]
                print ref_num, aln_num, seq_num
                out_file.write("%s %s %d %s %d %s" % (pdb, '\t', seq_num, '\t', ref_num, '\n'))
                