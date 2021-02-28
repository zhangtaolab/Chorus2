Chorus2 best practices
====================================

A guide pipeline for designing probes for Oligo-FISH in potato (*Solanum tuberosum*).
#############################################################################################

I. Requirements & Dependencies
------------------------------------------------

1.	Chorus2 software can be run on both Unix (MacOS/Linux) and Windows (via WSL) platform, here we recommend users to operate on Linux. Our best practices are also performed on an Ubuntu 18.04 system with 8 cores and 64 GB RAM.
2.	Chorus2 is written by python3, several packages are required: **Cython, numpy, pyfasta, primer3-py, pandas and pybigwig**. The GUI version of Chorus2 requires additional packages **matplotlib** and **PyQt5**.
3.	Besides, several bioinformatic tools are needed in Chorus2 pipeline. **BWA** is for genome-wide alignment. **Jellyfish** is for k-mer count.
4.	Genome assembly and whole genome sequencing (WGS) data are required for finding out oligos and filtering out repetitive sequences. WGS data with >5x depth are recommended for NGS filter process.
5.	Finally, we have already built a recipe on bioconda, thus Chorus2 can be installed easily via conda.


II. Preparation for Chorus2 pipeline
------------------------------------------------

1.	Download Anaconda distribution
*******************************************
Anaconda is an open-source distribution that is the easiest way to perform Python/R data science and machine learning on Linux. We use Anaconda3 to install Chorus2 directly.
Download the latest 64-bit Linux version of Anaconda3 distribution from the Anaconda website: https://www.anaconda.com/products/individual.

.. code-block:: bash

    $ cd /home/zhangtaolab/Software
    $ wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
    $ sh Anaconda3-2019.10-Linux-x86_64.sh
    $ export PATH=$HOME/anaconda3/bin:$PATH
    $ conda activate

After this process, software can be installed easily via conda command.

**Note**: The latest version of Anaconda3 was 2019.10 when we performed the best practices. If users want to use another version of Anaconda, please take a look from the official website.


2.	Install Chorus2 software
*******************************************
To install bioinformatic tools via conda, bioconda channels should be added.

.. code-block:: bash

    $ conda config --add channels defaults
    $ conda config --add channels bioconda

Then Chorus2 can be installed with conda command easily. We used conda to create an independent environment for Chorus2.

.. code-block:: bash

    $ conda create -n chorus Chorus2

After installing completed, check if Chorus2 runs well.
 
.. code-block:: bash

    $ conda activate chorus
    $ Chorus2 -h

This step will show the usage of Chorus2.


3.	Download reference genome file and whole genome sequencing files
************************************************************************
In this best practice, we designed oligos for potato genome. We used *Solanum tuberosum* species genome as input, which was sequenced by Potato Genome Sequence Consortium (PGSC). The reference genome file can be retrieved from the website: https://solgenomics.net/organism/Solanum_tuberosum/genome.

.. code-block:: bash

    $ cd /home/zhangtaolab/data/chorus2_project
    $ mkdir reference
    $ cd reference
    $ wget http://solanaceae.plantbiology.msu.edu/data/potato_dm_v404_all_pm_un.fasta.zip
    $ unzip potato_dm_v404_all_pm_un.fasta.zip

Whole genome sequencing data are used for repetitive sequences filtering after Chorus2 find out all available oligos. Here we used SRR5349606 data from NCBI Sequence Read Archive (SRA) database, which contains approximate 4x potato genome sequences.

.. code-block:: bash

    $ cd ..
    $ mkdir NGS
    $ cd NGS
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/006/SRR5349606/SRR5349606_1.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/006/SRR5349606/SRR5349606_2.fastq.gz

Here we used *wget* command to download SRR5349606 data, the links are copied from EBI database. Users can also use other tools to download WGS data, such as SRA Tookit and Aspera connect.


III. Perform Chorus2 pipeline for Oligo-FISH probes design
--------------------------------------------------------------------------------

1.	Using Chorus2 to design Oligo-FISH probes
******************************************************
Next, we can use Chorus2 software to design Oligo-FISH probes for potato.

.. code-block:: bash

    $ cd /home/zhangtaolab/data/chorus2_project
    $ mkdir analysis
    $ cd analysis
    $ Chorus2 -g /home/zhangtaolab/data/chorus2_project/reference/ potato_dm_v404_all_pm_un.fasta  -i /home/zhangtaolab/data/chorus2_project/reference/ potato_dm_v404_all_pm_un.fasta -s potato -t 16

Chorus2 program contains several parameters. '*-g*' requires a genome assembly file and '*-i*' requires an input file where oligos are designed. Both two are mandatory parameters. 
Other parameters are optional. '*-s*' requires a path where analysis results are saved. 
'*-p*' provides an examination of given 5' labeled R primer, which is used for PCR reaction. 
'*-t*' is used for speeding up analysis by add more threads. '*-l*' defines the fixed probe length. 
'*--homology*' defines the minimum homology between probe and target sequences. 
'*-d*' is for dTm check (hybrid *Tm* - hairpin *Tm*) to avoid duplex formation. 
'*--step'* defines the stride of k-mer searching in a sliding window.
Here our analysis results are stored in potato directory. File descriptions are as follow: 

*potato_dm_v404_all_pm_un.fasta.bed* is the probe file contained non-overlapped probes. 

*potato_dm_v404_all_pm_un.fasta_all.bed* is the probe file contained all probes. This file will be used for ChorusNGSfilter.

*potato_dm_v404_all_pm_un.fasta.len* is the length info of the given genome chromosomes. This file can be imported into ChorusPBGUI for probe selection.

*potato_dm_v404_all_pm_un.fasta_17mer.jf* is the binary file created by jellyfish count using 17-mer.

*potato_dm_v404_all_pm_un.fasta_tmp_probe.fa* contains all candidate probe sequences filtered by jellyfish.
*.bwt, .pac, .ann, .amb, .sa* files are bwa index files.

.. code-block:: bash

    $ cd /home/zhangtaolab/data/chorus2_project/analysis/potato
    $ head -5 potato_dm_v404_all_pm_un.fasta_all.bed
        chr01   97858   97902      ATTTTCCATGGACCTCATTAAGATTAGCTATTGAACCAGTTACCC
        chr01   103008  103052  ACAGCCAAATCGTCCCATATTCAAGGATAAACGACCCACGAATCA
        chr01   127095  127139  ATCTATATCTACTACACCAGAATATTCATACACAAATAAATTACT
        chr01   127101  127145  ATCTACTACACCAGAATATTCATACACAAATAAATTACTACTATT
        chr01   127815  127859  CTCAGATTTACCGAATTATTCCTGTGACAAAATATTACTTCCAGT

There are four columns in each row, first column is chromosome name, second is oligo start site, third is oligo end site, the last one is oligo probe sequence. Users can use excel or text editor to open this file.


2.	Use ChorusNGSfilter to filter repetitive sequences in oligos set
******************************************************************************
To further filter putative repetitive sequences, a kmer-based method can be performed to detect repeats by running ChorusNGSfilter. Here we used SRR5349606 data to facilitate filter process.

.. code-block:: bash

    $ ChorusNGSfilter -i /home/zhangtaolab/data/chorus2_project/NGS/SRR5349606_1.fastq.gz,/home/zhangtaolab/data/chorus2_project/NGS/SRR5349606_2.fastq.gz -z gz -g home/zhangtaolab/data/chorus2_project/reference/ potato_dm_v404_all_pm_un.fasta -t 16 -p potato_dm_v404_all_pm_un.fasta_all.bed -o potato_DM_v404_filtered.bed

ChorusNGSfilter requires WGS data as input (*-i*), genome file and designed probes file are also required(*-g* and *-p*). 
Length of k-mer used for counting k-mers can be adjusted by parameter '*-k*'.

After running NGS filtering, three files (*\*.jf, \*.bw, \*.bed*) will output to working directory:

*potato_DM_v404_filtered.bed.jf* is the binary file created by jellyfish count using given k-mer (Default is 17).

*potato_DM_v404_filtered.bed.bw* is a bigwig file contained all score information generated from NGS library.

*potato_DM_v404_filtered.bed* is the probe file contained all probes as well as k-mer score and strand. This file should be further selected by ChorusNGSselect.

.. code-block:: bash

    $ head -5 potato_dm_v404_all_pm_un.fasta_kmer.bed
        chr01   97858   97902   ATTTTCCATGGACCTCATTAAGATTAGCTATTGAACCAGTTACCC   324     +
        chr01   103008  103052  ACAGCCAAATCGTCCCATATTCAAGGATAAACGACCCACGAATCA   146     +
        chr01   127095  127139  ATCTATATCTACTACACCAGAATATTCATACACAAATAAATTACT   318     +
        chr01   127101  127145  ATCTACTACACCAGAATATTCATACACAAATAAATTACTACTATT   373     +
        chr01   127815  127859  CTCAGATTTACCGAATTATTCCTGTGACAAAATATTACTTCCAGT   281     +

There are six columns in each row, first four columns are the same as *Tpotato_dm_v404_all_pm_un.fasta_all.bed*. 
The fifth column is the k-mer score, last column is target strand of probes.


3.	Use ChorusNGSselect to select confident oligo probes
************************************************************************
Confident probes need to be further selected based on the k-mer of each oligo. ChorusNGSselect program is used for probes selection.

.. code-block:: bash

    $ ChorusNGSselect -i potato_DM_v404_filtered.bed -o potato_DM_v404_selected.bed
    $ head -5 potato_DM_v404_selected.bed
        chr01   97858   97902   ATTTTCCATGGACCTCATTAAGATTAGCTATTGAACCAGTTACCC   324     +
        chr01   103008  103052  TGATTCGTGGGTCGTTTATCCTTGAATATGGGACGATTTGGCTGT   146     -
        chr01   127095  127139  ATCTATATCTACTACACCAGAATATTCATACACAAATAAATTACT   318     +
        chr01   127815  127859  ACTGGAAGTAATATTTTGTCACAGGAATAATTCGGTAAATCTGAG   281     -
        chr01   133524  133568  ATTAATCAATAAAGGAAAAGCAAGGTTGGAATACGGTTTCATCCT   389     +

There are six columns in each row, which are the same as ChorusNGSfilter output. 
The final probes can be synthesized directly for oligo-FISH or imported into ChorusPBGUI for further selection.
