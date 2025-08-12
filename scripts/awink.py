from collections import defaultdict
import pickle as pkl
import csv
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gzip

def loadPickle(inFile):
    with open(inFile, 'rb') as file:
        data = pkl.load(file)
    return data

def writePickle(data, outFile):
    with open(outFile, 'wb') as file:
        pkl.dump(data, file)

def writeCSV(data, outFile):
    # Open the CSV file in write mode
    with open(outFile, mode='w', newline='') as file:
        # Create a CSV writer object
        writer = csv.writer(file)

        # Write each row of data to the CSV file
        for row in data:
            writer.writerow(row)

def readCSV(inFile):
    data = []
    with open(inFile, newline='') as file:
        reader = csv.reader(file)
        next(reader) # skip first row (header)
        for row in reader:
            data.append(row)
    return data

def get_read_to_unikmer_with_length_maps_summary(prefix):
    print('==============================')
    print(f"processing: started")
    read2UnikmerWithLengthMapSummary = dict()
    currentMap = loadPickle(prefix + ".pkl")
    for read in currentMap:
        read2UnikmerWithLengthMapSummary[read] = (currentMap[read]['length'], len(currentMap[read]['barcodes']))
    print(f"processing complete!!!")
    print('==============================')
    print(f'writing summary file!!!')
    writePickle(read2UnikmerWithLengthMapSummary, prefix + "_summary.pkl")
    print(f"writing complete!!!")

def get_read_to_unikmer_with_length_maps_partial_summary(prefix, b):
    print('==============================')
    print(f"processing: {b}")
    read2UnikmerWithLengthMapSummary = dict()
    currentMap = loadPickle(prefix + "_" + str(b) + ".pkl")
    for read in currentMap:
        read2UnikmerWithLengthMapSummary[read] = (currentMap[read]['length'], len(currentMap[read]['barcodes']))
    print(f"=== testing === {b} ===")
    tmp_read = list(read2UnikmerWithLengthMapSummary.keys())[0]
    print(f"readID: {tmp_read} -> {read2UnikmerWithLengthMapSummary[tmp_read]}")
    print(f"processing complete!!!")
    print('==============================')
    print(f'writing summary file!!!')
    writePickle(read2UnikmerWithLengthMapSummary, prefix + "_" + str(b) + "_summary.pkl")
    print(f"writing complete!!!")

def join_read_to_unikmer_with_length_maps_summary(prefix, bin):
    read2UnikmerWithLengthMapSummary = dict()
    for b in range(1,bin+1):
        print('==============================')
        print(f"processing: {b}")
        currentMap = loadPickle(prefix + "_" + str(b) + "_summary.pkl")
        read2UnikmerWithLengthMapSummary = {**read2UnikmerWithLengthMapSummary, **currentMap}
        print(f"processing complete!!!")
    print('==============================')
    print(f'writing summary file!!!')
    writePickle(read2UnikmerWithLengthMapSummary, prefix + "_all_summary.pkl")
    print(f"writing complete!!!")

def get_reads_dict(input_fastq):
    print("started reading...")
    reads_dict = dict()
    reads = ""
    if input_fastq.endswith(".gz"):
        # Open the FASTQ.gz file for reading
        with gzip.open(input_fastq, "rt") as handle:
            # Iterate over each record (sequence) in the FASTQ/FASTA file
            if input_fastq.endswith(".fastq.gz"):
                reads = SeqIO.parse(handle, "fastq")
            elif input_fastq.endswith(".fasta.gz"):
                reads = SeqIO.parse(handle, "fasta")
            for record in reads:
                # Access each sequence using the 'seq' attribute of the record object
                reads_dict[record.id] = record.seq
    else:
        # Open the FASTQ file for reading    
        with open(input_fastq, "r") as handle:
            # Iterate over each record (sequence) in the FASTQ/FASTA file
            if input_fastq.endswith(".fastq"):
                reads = SeqIO.parse(handle, "fastq")
            elif input_fastq.endswith(".fasta"):
                reads = SeqIO.parse(handle, "fasta")
            for record in reads:
                # Access each sequence using the 'seq' attribute of the record object
                reads_dict[record.id] = record.seq
    print("reading completed: number of reads ->", len(reads_dict))
    return reads_dict

def get_unikmer_distribution_per_read(readBarcodeWithLengthMapFile, unikmerDistributionFile, outCSVFile):
    print(f"1. === loading read barcode with length summary file...")
    readBarcodeWithLengthMap = loadPickle(readBarcodeWithLengthMapFile)

    print(f"2. === processing unikmer distribution over reads...")
    unikmerDistribution = defaultdict(dict)
    for read in readBarcodeWithLengthMap:
        unikmerDistribution[readBarcodeWithLengthMap[read][1]][read] = readBarcodeWithLengthMap[read][0]
    sortedUnikmerDistribution = dict(sorted(unikmerDistribution.items()))

    print(f"3. === writing sorted unikmer distribution map...")
    writePickle(sortedUnikmerDistribution, unikmerDistributionFile)

    unikmerDistributionOverReads = [["# of unikmers", "# of reads", "sum of bases (bp)"]]
    for unikmerCount in sortedUnikmerDistribution:
        baseSum = 0
        for read in sortedUnikmerDistribution[unikmerCount]:
            baseSum += sortedUnikmerDistribution[unikmerCount][read]
        unikmerDistributionOverReads.append([unikmerCount, len(sortedUnikmerDistribution[unikmerCount]), baseSum])
    
    print(f"4. === writing unikmer distribution to csv file...")
    writeCSV(unikmerDistributionOverReads, outCSVFile)
    
    print(f"5. === done!!!")


# codes for sampling following the unikmer distribution prioritizing long reads
def get_mean_std_for_unikmer_distribution_per_read(inputCSVFile):
    no_of_unikmers, no_of_reads = [], []
    with open(inputCSVFile, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)

        for row in reader:
            no_of_unikmers.append(int(row[0]))
            no_of_reads.append(int(row[1]))
    
    start_index = 50
    print(f"start index: {start_index}")
    truncated_no_of_reads = no_of_reads[start_index:]
    max_index = start_index + truncated_no_of_reads.index(max(truncated_no_of_reads))
    print(f"max value: {max(truncated_no_of_reads)}")
    print(f"max index: {max_index}")
    end_index = start_index + 2*(max_index - start_index)
    print(f"end index: {end_index}")

    # calculating the mean and the standard deviation of the distribution

    F = np.array(no_of_reads[start_index:])
    X = np.array(no_of_unikmers[start_index:])

    mean = np.sum(X*F)/np.sum(F)
    variance = np.sum(F*((X-mean)*(X-mean)))/np.sum(F)
    stdev = np.sqrt(variance)

    print(f"mean: {mean} std: {stdev}")
    
    # calculating the lower and upper bounds for extracting unikmers
    
    lower = max(int(np.floor(mean - 3*stdev)), start_index)
    upper = int(np.ceil(mean + 3*stdev))

    print(f"lower: {lower}")
    print(f"upper: {upper}")

    return lower, upper, no_of_unikmers[0], no_of_unikmers[-1]

# adaptive binning based on full coverage and genome size
# for every 1000x full coverage of a 0.5-5 MB genome: 10x per bin
# for every 1000x full coverage of a 5-50 MB genome: 1x per bin
# for every 1000x full coverage of a 50-500 MB genome: 0.1x per bin
# for every 1000x full coverage of a 500 MB or larger genome: 0.01x per bin
def adaptive_read_binning_based_on_unikmer_distribution(readBarcodeWithLengthMapFile, inputCSVFile, binnedReadsOnUnikmerCountFile, coverage, genomeSize, outputCSVFile):
    lower, upper, freq_min, freq_max = get_mean_std_for_unikmer_distribution_per_read(inputCSVFile)
    print(f"min: {freq_min}, max: {freq_max}")
    print(f"lower: {lower}, upper: {upper}")

    # making it a multiple of 100
    adjusted_lower = int(np.round(lower/100))*100
    # making it a multiple of 500
    adjusted_upper = int(np.round(upper/500))*500

    print(f"adjusted lower: {adjusted_lower}, adjusted upper: {adjusted_upper}")

    print(f"1. === loading read barcode with length summary file...")
    readBarcodeWithLengthMap = loadPickle(readBarcodeWithLengthMapFile)

    print(f"2. === binning reads based on unikmer counts...")
    
    # determining minimum coverage per bin
    
    coveragePerBin = -1.0
    if 0 <= genomeSize < 500000:
        print("too small genome for binning reads")
        return 
    elif 500000 <= genomeSize < 5000000:
        coveragePerBin = coverage / 100
    elif 5000000 <= genomeSize < 50000000:
        coveragePerBin = coverage / 1000
    elif 50000000 <= genomeSize < 500000000:
        coveragePerBin = coverage / 10000
    elif 500000000 >= genomeSize:
        coveragePerBin = coverage / 100000

    minimumBaseSumPerBin = int(genomeSize * coveragePerBin)
    unikmerDistributionCSV = readCSV(inputCSVFile)
    baseSumPerUnikmer = dict()
    unikmerCountList = list()
    for row in unikmerDistributionCSV:
        baseSumPerUnikmer[int(row[0])] = int(row[2])
        unikmerCountList.append(int(row[0]))

    binnedReadsOnUnikmerCount = dict()
    unikmer2BinMap = dict()

    # create bins
    # always put reads with unikmers 0-5 in separate bins
    currentIndex = 0
    for i in range(6):
        if i in baseSumPerUnikmer:
            binnedReadsOnUnikmerCount[i] = dict()
            unikmer2BinMap[i] = i
            currentIndex += 1
    
    currentWindowStart = unikmerCountList[currentIndex]
    currentWindowBaseSum = 0
    currentUnikmerIndex = currentIndex
    # for currentUnikmer in unikmerCountList[currentIndex:]:
    while unikmerCountList[currentUnikmerIndex] <= adjusted_upper:
        currentWindowBaseSum += baseSumPerUnikmer[unikmerCountList[currentUnikmerIndex]]
        unikmer2BinMap[unikmerCountList[currentUnikmerIndex]] = currentWindowStart
        # print(f"{currentUnikmer} : {currentWindowStart}")
        if currentWindowBaseSum >= minimumBaseSumPerBin:
            binnedReadsOnUnikmerCount[currentWindowStart] = dict()
            currentWindowBaseSum = 0
            currentUnikmerIndex += 1
            if currentUnikmerIndex < len(unikmerCountList):
                currentWindowStart = currentUnikmerIndex
        else:
            currentUnikmerIndex += 1
    
    if currentWindowStart not in binnedReadsOnUnikmerCount:
        binnedReadsOnUnikmerCount[currentWindowStart] = dict()

    binnedReadsOnUnikmerCount[adjusted_upper+1] = dict()

    for read in readBarcodeWithLengthMap:
        length, unikmerCount = readBarcodeWithLengthMap[read][-2], readBarcodeWithLengthMap[read][-1]
        if 0 <= unikmerCount <= 5:
            binnedReadsOnUnikmerCount[unikmerCount][read] = length
        elif unikmerCount > adjusted_upper:
            binnedReadsOnUnikmerCount[adjusted_upper+1][read] = length
        else:
            binnedReadsOnUnikmerCount[unikmer2BinMap[unikmerCount]][read] = length

    for slab in binnedReadsOnUnikmerCount:
        temp_dict = dict(sorted(binnedReadsOnUnikmerCount[slab].items(), key=lambda x:x[1], reverse=True))
        binnedReadsOnUnikmerCount[slab] = temp_dict

    print(f"3. === writing read bins based on unikmer counts to pkl file...")
    writePickle(binnedReadsOnUnikmerCount, binnedReadsOnUnikmerCountFile)

    print(f"4. === writing read bins based on unikmer counts to csv file...")
    unikmerDistributionOverReads = [["unikmer bin", "# of reads", "sum of bases (bp)"]]
    for slab in binnedReadsOnUnikmerCount:
        baseSum = 0
        for read in binnedReadsOnUnikmerCount[slab]:
            baseSum += binnedReadsOnUnikmerCount[slab][read]
        unikmerDistributionOverReads.append([slab, len(binnedReadsOnUnikmerCount[slab]), baseSum])

    writeCSV(unikmerDistributionOverReads, outputCSVFile)
    
    print(f"5. === done!!!")


def select_reads_based_on_unikmer_distribution(binnedReadsOnUnikmerCountFile, targetCoverages, coverage, genomeSize, selectedReadsFile, binnedReadsOnUnikmerCountCSVFile):
    binnedReadsOnUnikmerCount = loadPickle(binnedReadsOnUnikmerCountFile)
    selectedReadIDs = dict()
    
    binnedReadsOnUnikmerCountCSV = readCSV(binnedReadsOnUnikmerCountCSVFile)
    last_slab, no_of_reads_last_slab, sum_of_bases_last_slab = int(binnedReadsOnUnikmerCountCSV[-1][0]), int(binnedReadsOnUnikmerCountCSV[-1][1]), int(binnedReadsOnUnikmerCountCSV[-1][2])
    initial_coverage = sum_of_bases_last_slab / genome_size
    samplingRatios = [(target - initial_coverage)/coverage for target in targetCoverages]

    binnedReadsOnUnikmerCountSamplingBaseSums = dict()
    for slab in binnedReadsOnUnikmerCount:
        baseSum = 0
        for read in binnedReadsOnUnikmerCount[slab]:
            baseSum += binnedReadsOnUnikmerCount[slab][read]
        binnedReadsOnUnikmerCountSamplingBaseSums[slab] = [int(ratio*baseSum) for ratio in samplingRatios]
    binnedReadsOnUnikmerCountSamplingBaseSums[last_slab] = [0 for ratio in samplingRatios]

    for i in range(len(targetCoverages)):
        selectedReadIDs[targetCoverages[i]] = set(binnedReadsOnUnikmerCount[last_slab].keys())
        # print(len(selectedReadIDs[targetCoverages[i]]))
        for slab in binnedReadsOnUnikmerCount:
            currentSum = 0
            currentReads = list(binnedReadsOnUnikmerCount[slab].keys())
            # print(slab, currentReads[:3])
            j = 0
            while currentSum < binnedReadsOnUnikmerCountSamplingBaseSums[slab][i] and j < len(currentReads):
                selectedReadIDs[targetCoverages[i]].add(currentReads[j])
                currentSum += binnedReadsOnUnikmerCount[slab][currentReads[j]]
                j += 1

    writePickle(selectedReadIDs, selectedReadsFile)        


def extract_selected_binned_reads(selectedReadsFile, readFile, outFile, depths):
    selectedReads = loadPickle(selectedReadsFile)
    readsDict = get_reads_dict(readFile)
    prevReads = set()
    for i in range(len(depths)):
        print("==============================")
        print(f"current depth: {depths[i]}x")
        readIDs = list(selectedReads[depths[i]] - prevReads)
        print(f"number of reads: {len(readIDs)}")
        seqs = []
        for readID in readIDs:
            seqs.append(SeqRecord(readsDict[readID], readID))
        print(f"writing current reads...")
        SeqIO.write(seqs, outFile + str(depths[i]) + "x.fasta", "fasta")
        prevReads = selectedReads[depths[i]]
    print("==============================")
    print("writing extracted reads: completed!!!")


def awink_executor():
    k = "value of k for k-mer; in this study k=21"
    coverage = "coverage depth of the ultra-deep sequencing data"

    path = "path to parent directory"
    genome_size = "provide the estimated genome size"

    coverage_depths = "provide a list of target coverage depths as input to AWinK" # example: [20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]

    readFile = path + "/path to ultra-deep sequencing read file (fastq/fasta)"

    readBarcodeWithLengthMapFile = path + "/barcode_"+str(k)+"mers/read2UnikmerWithReadLengthMap_summary.pkl"
    readBarcodeWithLengthMapFilePath = path + "/barcode_"+str(k)+"mers/read2UnikmerWithReadLengthMap"
    
    unikmerDistributionFile = path + "/barcode_"+str(k)+"mers/unikmerDistributionMap.pkl"
    unikmerDistributionCSVFile = path + "/barcode_"+str(k)+"mers/unikmerDistribution.csv"
    binnedReadsOnUnikmerCountFile = path + "/barcode_"+str(k)+"mers/binnedReadsOnUnikmerCount_awink.pkl"
    binnedReadsOnUnikmerCountCSVFile = path + "/barcode_"+str(k)+"mers/binnedReadsOnUnikmerCount_awink.csv"
    selectedBinnedReadIDsFile = path + "/barcode_"+str(k)+"mers/reads/selected_awink" + "_" + str(coverage_depths[-1])+"x.pkl"
    extractedBinnedReadsFilePath = path + "/barcode_"+str(k)+"mers/reads/extracted_awink_"


    # function calls for extracting reads based on adaptive binning from barcoded reads
    print(f"step 1: get read barcode summary")
    get_read_to_unikmer_with_length_maps_summary(readBarcodeWithLengthMapFilePath)
    print(f"step 2: generate unikmer distribution per read barcode")
    get_unikmer_distribution_per_read(readBarcodeWithLengthMapFile, unikmerDistributionFile, unikmerDistributionCSVFile)
    print(f"step 3: perform read binning based on adaptive window")
    adaptive_read_binning_based_on_unikmer_distribution(readBarcodeWithLengthMapFile, unikmerDistributionCSVFile, binnedReadsOnUnikmerCountFile, coverage, genome_size, binnedReadsOnUnikmerCountCSVFile)
    print(f"step 4: select reads from bins")
    select_reads_based_on_unikmer_distribution(binnedReadsOnUnikmerCountFile, coverage_depths, coverage, genome_size, selectedBinnedReadIDsFile, binnedReadsOnUnikmerCountCSVFile)
    print(f"step 5: extract selected reads")
    extract_selected_binned_reads(selectedBinnedReadIDsFile, readFile, extractedBinnedReadsFilePath, coverage_depths)


if __name__ == "__main__":
    awink_executor()