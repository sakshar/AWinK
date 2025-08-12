from Bio import SeqIO
import pickle as pkl
import os

def loadPickle(inFile):
    with open(inFile, 'rb') as file:
        data = pkl.load(file)
    return data

def writePickle(data, outFile):
    with open(outFile, 'wb') as file:
        pkl.dump(data, file)

def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def get_unikmer_map(unikmer_file):
    unikmer_map = dict()
    counter = 0
    with open(unikmer_file, "r") as file:
        for line in file:
            unikmer_map[line.strip()] = counter
            counter += 1
    return unikmer_map

def barcode_reads_with_length(unikmer_map, read_file, k):
    read_to_unikmer_with_readLength_map = dict()
    with open(read_file, "r") as handle:
        # Iterate over each record (sequence) in the FASTQ/FASTA file
        reads = ""
        if read_file.endswith(".fastq"):
            reads = SeqIO.parse(handle, "fastq")
        elif read_file.endswith(".fasta"):
            reads = SeqIO.parse(handle, "fasta")
        counter = 0
        for record in reads:
            read = record.seq
            id = record.id
            read_to_unikmer_with_readLength_map[id] = {'length': len(read), 'barcodes': list()}
            for i in range(len(read) - k + 1):
                kmer = read[i:i + k]
                kmer_rev = kmer.reverse_complement()
                if kmer in unikmer_map:
                    read_to_unikmer_with_readLength_map[id]['barcodes'].append(unikmer_map[kmer])
                elif kmer_rev in unikmer_map:
                    read_to_unikmer_with_readLength_map[id]['barcodes'].append(unikmer_map[kmer_rev])
            counter += 1
            if counter % 50000 == 0:
                print(f"processed number of reads: {counter}")
    return read_to_unikmer_with_readLength_map


def mapper_executor():
    k = "value of k for k-mer; in this study k=21"
    depth = "coverage depth of the ultra-deep sequencing data"

    dir = "path to parent directory"
    unikmerTxtFile = dir + "path to $prefix.unikmers inside the jellyfish output directory"
    readFile = dir + "/path to ultra-deep sequencing read file (fastq/fasta)"

    outDir = dir + "/barcode_"+str(k)+"mers"
    create_directory(outDir)

    unikmerMapFile = outDir + "/unikmerMap.pkl"
    readBarcodeWithLengthMapFile = outDir + "/read2UnikmerWithReadLengthMap.pkl"

    print("====== initializing ======")
    print(f"input HiFi reads: {readFile}")
    print(f"output map: {readBarcodeWithLengthMapFile}")
    print("==========================")
    print(f"step 1: loading unikmers")
    unikmerMap = get_unikmer_map(unikmerTxtFile)
    writePickle(unikmerMap, unikmerMapFile)
    print(f"step 2: barcoding reads with unikmers and lengths")
    read2unikmerLengthMap = barcode_reads_with_length(unikmerMap, readFile, k)
    print(f"step 3: writing read barcodes with read lengths")
    writePickle(read2unikmerLengthMap, readBarcodeWithLengthMapFile)
    print("======== done!!! =========")


if __name__ == "__main__":
    mapper_executor()
