# AWinK

This is a stand-alone repository for executing AWinK in your local machine. AWinK is developed as a part of a project on improving de novo genome assembly in presence of ultra-deep sequencing data.

# Pre-requisites
You will need conda installed in your machine to proceed. If you don't have conda installed, please follow [conda installation](https://conda.io/docs/user-guide/install/). If you already had it installed, please make sure it is updated.

Create a conda environment first:

`conda create -n awink python=3.8 numpy xlsxwriter biopython`

This will create a conda environment named rambler with the afforementioned packages. Then activate the environment with:

`conda activate awink`

Finally, install the following packages inside the newly created environment using the commands below:

<pre>
conda install -c bioconda jellyfish
conda install -c bioconda minimap2
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda pysam
conda install -c bioconda hifiasm
conda install -c bioconda seqkit
</pre>

# Scripts

`extractUnikmers.sh` contains the required commands to perform step A (Extract unikmers)
`read2UnikmerMap.py` contains the required functions to perform step B (Barcode reads)  
`awink.py` contains the required fucntions to execute AWinK's steps C and D (Assign reads to bins and Select reads)
`assembly.sh` contains the required commands to assemble the AWinK selected reads using hifiasm (step E: Assemble reads)

# Execution

Open the terminal and follow the steps below to run RAmbler from scratch:
<pre>
  cd ~
  git clone https://github.com/sakshar/AWinK.git
  cd AWinK/scripts
  bash extractUnikmers.sh prefix k
  python read2UnikmerMap.py
  python awink.py
  bash assembly.sh
</pre>

Executing the above commands one-by-one will generate the final assembly inside the user-specified hifiasm output directory with the filename `$PREFIX.asm.fasta`
