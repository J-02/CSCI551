from Bio import SeqIO
from Bio.Seq import Seq
def getSeq(file):
    record = SeqIO.read(file, "fasta")
    seq = Seq(str(record.seq)+ "$")
    return seq
def fmindex(string):
