from Bio import SeqIO
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Custom fastq trimmer for Illumina-Sanger')
    parser.add_argument('-i', help='Input fastq file', metavar='Str', type=str)
    parser.add_argument('-he', help='Number of nucleotides cut from the beginning og the read', metavar = 'Int', type = int, default=0)
    parser.add_argument('-c', help='Number of nucleotides cut from the end of the sequence', metavar = 'Int', type = int, default=0)
    parser.add_argument('-swl', help='Sliding window length', metavar = 'Int', type = int, default = 0)
    parser.add_argument('-swq', help='Sliding window phred33 quality baseline', metavar = 'Int', type = int, default = 0)
    parser.add_argument('-o', help='Output fastq file', metavar='Str', type=str)
    args = parser.parse_args()
    inp = args.i
    outp = args.o
    head = args.he
    crop = args.c
    length = args.swl
    qual = args.swq

def HEADCROP(read, drop_value):
#    new_seq, new_an = read.seq[drop_value - 1:], {'phred_quality': read.letter_annotations['phred_quality'][drop_value - 1:]}
#    new_read = SeqIO.SeqRecord(seq=new_seq, letter_annotations=new_an, id=read.id, features=read.features)
#hre and hereafter we will use a simpler variant
    new_read = read[drop_value - 1:]
    return new_read

def CROP(read, drop_value):
    new_read = read[:len(read) - drop_value]
    return new_read

def SLIDING_WINDOW(read, length, baseline):
    start = 0
    for i in range(start, len(read) + 1 - length):
        start = i
        if sum(read.letter_annotations['phred_quality'][start:start + length])/length < baseline:
            new_read = CROP(read, start)
            return new_read

def trimming(file, output):
    list = []
    for record in SeqIO.parse(file, 'fastq'):
        a = CROP(record, head)
        if not (a is None):
            b = HEADCROP(a, crop)
            if (not b is None):
                c = SLIDING_WINDOW(b, length, qual)
                if not (c is None) and len(c.seq) > 10:
                    list.append(c)
    SeqIO.write(list, output, 'fastq')

trimming(inp, outp)