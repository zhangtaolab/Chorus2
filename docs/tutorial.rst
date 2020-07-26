Tutorial
=========

Use Chorus2 to design oligo probes for plant genome
------------------------------------------------------------

In this tutorial, we will build oligo probe set for Arabidopsis genome.

Install Chorus2
******************************

See install tutorial here_

.. _here: install.html

Run Chorus2
******************************

Run Chorus2 with Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Download Reference Genome file**

.. code-block:: bash

    $ wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

    $ docker run -v $PWD:/home/chorus -e CHORUS_USER=$USER -e CHORUS_UID=$UID \ 
      forrestzhang/docker-chorus -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 12

Please wait unit all precess done. There are some logs:

.. code-block:: log

    forrest /home/chorus
    use local user:  forrest
    Adding group 'forrest' (GID 1000) ...
    Done.
    Adding user 'forrest' ...
    Adding new user 'forrest' (1000) with group 'forrest' ...
    Creating home directory '/home/forrest' ...
    Copying files from '/etc/skel' ...
    /home/chorus exists
    2.2.3
    ########################################
    bwa version: /opt/software/bwa/bwa 0.7.12-r1044
    jellyfish version: /opt/software/jellyfish/bin/jellyfish 2.2.3
    genome file: TAIR10_chr_all.fas
    input file: TAIR10_chr_all.fas
    5' labeled R primer:
    result output folder: /home/chorus/probes
    threads number: 12
    homology: 75
    dtm: 10
    ########################################
    ...
    ...
    14300000 / 14326857
    14310000 / 14326857
    14320000 / 14326857
    Job finshed!!


When process done:

.. code-block:: bash

    $ ls -lt probes/
    total 1741428
    -rw-r--r-- 1 root root  280927981 Aug 24 17:44 TAIR10_chr_all.fas_all.bed
    -rw-r--r-- 1 root root   62050561 Aug 24 17:44 TAIR10_chr_all.fas.bed
    -rw-r--r-- 1 root root         94 Aug 24 17:30 TAIR10_chr_all.fas.len
    -rw-r--r-- 1 root root 1031512169 Aug 24 17:22 TAIR10_chr_all.fas_tmp_probe.fa
    -rw-r--r-- 1 root root   59833928 Aug 24 17:19 TAIR10_chr_all.fas.sa
    -rw-r--r-- 1 root root       7535 Aug 24 17:18 TAIR10_chr_all.fas.amb
    -rw-r--r-- 1 root root        682 Aug 24 17:18 TAIR10_chr_all.fas.ann
    -rw-r--r-- 1 root root   29916939 Aug 24 17:18 TAIR10_chr_all.fas.pac
    -rw-r--r-- 1 root root  119667836 Aug 24 17:18 TAIR10_chr_all.fas.bwt
    -rw-r--r-- 1 root root  121183059 Aug 24 17:17 TAIR10_chr_all.fas
    -rw-r--r-- 1 root root   78102510 Aug 24 17:17 TAIR10_chr_all.fas_17mer.jf

*TAIR10_chr_all.fas.bed* is the probe file contained non-overlapped probes.

*TAIR10_chr_all.fas_all.bed* is the probe file contained all probes. This file can be used for ChorusNGSfilter.

*TAIR10_chr_all.fas.len* is the length info of the given genome chromosomes. This file can be imported into ChorusPBGUI for probe selection.

*TAIR10_chr_all.fas_17mer.jf* is the binary file created by jellyfish count using 17-mer.

*TAIR10_chr_all.fas_tmp_probe.fa* contains all candidate probe sequences filtered by jellyfish.

*.bwt, .pac, .ann, .amb, .sa* files are bwa index files.

.. code-block:: log

    $ more probes/TAIR10_chr_all.fas.bed
    1      	52     	96     	TCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTA
    1      	211    	255    	TTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTC
    1      	346    	390    	CCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACT
    1      	426    	470    	TTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAAT
    1      	496    	540    	TCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGG
    1      	551    	595    	TAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAG

There are four columns in each row, first column is chromosome name, second is oligo start site, third is oligo end site, the last one is oligo probe sequence. You can use excel or text editor to open this file.


Run Chorus2 in terminal
^^^^^^^^^^^^^^^^^^^^^^^^

**Make a project folder**

.. code-block:: bash

    $ cd ~
    $ mkdir sampleproject
    $ cd sampleproject

**Download reference genome**

.. code-block:: bash

    $ wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas


**Test chorus2 software**

.. code-block:: bash

    $ Chorus2 -h
    usage: Chorus2 [-h] [--version] [-j JELLYFISH] [-b BWA] -g GENOME -i INPUT
                [-s SAVED] [-p PRIMER] [-t THREADS] [-l LENGTH]
                [--homology HOMOLOGY] [-d DTM] [--skipdtm SKIPDTM]
                [--step STEP] [--docker DOCKER] [--ploidy PLOIDY]

    Chorus2 Software for Oligo FISH probe design

    optional arguments:
    -h, --help            show this help message and exit
    --version             show program\'s version number and exit
    -j JELLYFISH, --jellyfish JELLYFISH
                          The path where Jellyfish software installed
    -b BWA, --bwa BWA     The path where BWA software installed
    -g GENOME, --genome GENOME
                          Fasta format genome file, should include all sequences
                          from genome
    -i INPUT, --input INPUT
                          Fasta format input file, can be whole genome, a
                          chromosome or one region from genome
    -s SAVED, --save SAVED
                          The output folder for saving results
    -p PRIMER, --primer PRIMER
                          A specific 5\' labeled R primer for PCR reaction. For
                          example: CGTGGTCGCGTCTCA. (Default is none)
    -t THREADS, --threads THREADS
                          Number of threads or CPUs to use. (Default: 1)
    -l LENGTH, --length LENGTH
                          The probe length. (Default: 45)
    --homology HOMOLOGY   The maximum homology(%) between target sequence and
                          probe, range from 50 to 100. (Default: 75)
    -d DTM, --dtm DTM     The minimum value of dTm (hybrid Tm - hairpin Tm),
                          range from 0 to 37. (Default: 10)
    --skipdtm SKIPDTM     skip calculate dtm, for oligo longer than 50.
    --step STEP           The step length for k-mer searching in a sliding
                          window, step length>=1. (Default: 5)
    --docker DOCKER       Only used in Docker version of Chorus
    --ploidy PLOIDY       The ploidy of the given genome (test version).
                          (Default: 2)

    Example:
    Chorus2 -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 4 \
            -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa -s sample

**Run chorus2 software**

.. code-block:: bash

    $ Chorus2 -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 12

When job finish, the oligo probes will output to 'probes' folder (Default, can be changed using -s)

.. code-block:: bash

    $ cd sample
    $ ls -lt * 

        total 1741428
        -rw-r--r-- 1 root root  280927981 Aug 24 17:44 TAIR10_chr_all.fas_all.bed
        -rw-r--r-- 1 root root   62050561 Aug 24 17:44 TAIR10_chr_all.fas.bed
        -rw-r--r-- 1 root root         94 Aug 24 17:30 TAIR10_chr_all.fas.len
        -rw-r--r-- 1 root root 1031512169 Aug 24 17:22 TAIR10_chr_all.fas_tmp_probe.fa
        -rw-r--r-- 1 root root   59833928 Aug 24 17:19 TAIR10_chr_all.fas.sa
        -rw-r--r-- 1 root root       7535 Aug 24 17:18 TAIR10_chr_all.fas.amb
        -rw-r--r-- 1 root root        682 Aug 24 17:18 TAIR10_chr_all.fas.ann
        -rw-r--r-- 1 root root   29916939 Aug 24 17:18 TAIR10_chr_all.fas.pac
        -rw-r--r-- 1 root root  119667836 Aug 24 17:18 TAIR10_chr_all.fas.bwt
        -rw-r--r-- 1 root root  121183059 Aug 24 17:17 TAIR10_chr_all.fas
        -rw-r--r-- 1 root root   78102510 Aug 24 17:17 TAIR10_chr_all.fas_17mer.jf

*TAIR10_chr_all.fas.bed* is the probe file contained non-overlapped probes.

*TAIR10_chr_all.fas_all.bed* is the probe file contained all probes. This file can be used for ChorusNGSfilter.

*TAIR10_chr_all.fas.len* is the length info of the given genome chromosomes. This file can be imported into ChorusPBGUI for probe selection.

*TAIR10_chr_all.fas_17mer.jf* is the binary file created by jellyfish count using 17-mer.

*TAIR10_chr_all.fas_tmp_probe.fa* contains all candidate probe sequences filtered by jellyfish.

*.bwt, .pac, .ann, .amb, .sa* files are bwa index files.

.. code-block:: log

    $ more probes/TAIR10_chr_all.fas.bed
    1      	52     	96     	TCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTA
    1      	211    	255    	TTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTC
    1      	346    	390    	CCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACT
    1      	426    	470    	TTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAAT
    1      	496    	540    	TCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGG
    1      	551    	595    	TAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAG

There are four columns in each row, first column is chromosome name, second is oligo start site, third is oligo end site, the last one is oligo probe sequence. You can use excel or text editor to open this file.

**Further filter using ChorusNGSfilter**

To further filter putative repetitive sequences, a kmer-based method can be performed to detect repeats by running ChorusNGSfilter.
Before running ChorusNGSfilter, a set of whole-genome shotgun sequencing data is required. Here we download the shotgun reads of Arabidopsis with the accession number SRR5658649.

.. code-block:: bash

    $ wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR565/009/SRR5658649/SRR5658649_1.fastq.gz
    $ wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR565/009/SRR5658649/SRR5658649_2.fastq.gz

    $ ChorusNGSfilter -i SRR5658649_1.fq.gz,SRR5658649_2.fq.gz -z gz \
                      -g TAIR10_chr_all.fas -t 12 \
                      -p probes/TAIR10_chr_all.fas_all.bed -o probes/TAIR10_chr_all_SRR5658649.bed

After running NGS filtering, three files (\*.jf, \*.bw, \*.bed) will output to working directory.

*TAIR10_chr_all_SRR5658649.bed.jf* is the binary file created by jellyfish count using given k-mer (Default is 17).

*TAIR10_chr_all_SRR5658649.bed.bw* is a bigwig file contained all score infomation generated from NGS library.

*TAIR10_chr_all_SRR5658649.bed* is the the probe file contained all probes as well as k-mer score and strand. This file should be further selected by ChorusNGSselect.

.. code-block:: log

    $ more probes/TAIR10_chr_all_SRR5658649.bed
    1	12	56	AAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCT	455128	+
    1	18	62	TAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCT	346	    +
    1	24	68	CTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAAT	343	    +
    1	36	80	ATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATG	319	    +
    1	42	86	AATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCC	315	    +
    1	48	92	TAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATA	294	    +

There are six columns in each row, first four columns are the same as TAIR10_chr_all.fas_all.bed. The fifth column is the k-mer score, last column is target strand of probes.

**Automatic probe selection using ChorusNGSselect**

Probes should be filtered by kmer score, the process can be done by ChorusNGSselect.

.. code-block:: bash

    $ ChorusNGSselect -i probes/TAIR10_chr_all_SRR5658649.bed \
                      -o probes/TAIR10_chr_all_SRR5658649_filter.bed

ChorusNGSselect will generate a final filtered probe file, it looks this:

.. code-block:: log

    $ more probes/TAIR10_chr_all_SRR5658649_filter.bed
    1    36      80      ATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATG   319     +
    1    66      110     CGGGTTTAGGGAATTAGGTATTTAGGGATTCATGGATGTAGGATT   221     -
    1    215     259     AGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTCACTT   293     +
    1    245     289     ATAACAAATGAAGATAAACCATCCATAGCTAAGTGAAGGAAAGAA   291     -
    1    347     391     CTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTG   237     +
    1    425     469     TTTTCCCACAAAAGATAAATTGATAGAACAAACATTTCCACAAAG   360     -

The final probes can be synthesized directly for oligo-FISH or imported into ChorusPBGUI for further selection.


Run Chorus2 with GUI
^^^^^^^^^^^^^^^^^^^^^^^^

**Make a project folder**

.. code-block:: bash

    $ cd ~
    $ mkdir sampleproject
    $ cd sampleproject

**Download reference genome**

.. code-block:: bash

    $ wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

**Run ChorusGUI**

.. code-block:: bash

    $ ChorusGUI

Set your own parameters and click Run to start the design process.

When job finish, the oligo probes will output to Sample Folder where you set.

**Further filter using ChorusNGSfilter**

The same process as "**Run Chorus2 in terminal**"

**Automatic probe selection using ChorusNGSselect**

The same process as "**Run Chorus2 in terminal**"

**Run ChorusPBGUI**

After filtering the probes, users can select suitable number of probes in specific regions 
for their FISH experiments using ChorusPBGUI easily.

.. code-block:: bash

    $ ChorusPBGUI



Use ChorusHomo to design oligo probes for close related species
-----------------------------------------------------------------

In this tutorial, we will design oligo probes for potato genome and detect conserved oligos 
in its close related species, tomato.

Run ChorusHomo
******************************

**Make a project folder**

.. code-block:: bash

    $ cd ~
    $ mkdir sampleproject
    $ cd sampleproject


**Download Reference Genome file**

.. code-block:: bash

    $ wget http://solanaceae.plantbiology.msu.edu/data/potato_dm_v404_all_pm_un.fasta.zip
    $ unzip potato_dm_v404_all_pm_un.fasta.zip
    $ wget ftp://ftp.solgenomics.net/tomato_genome/assembly/build_3.00/S_lycopersicum_chromosomes.3.00.fa


**Run Chorus2 to design probes in source species**

.. code-block:: bash

    $ Chorus2 -i potato_dm_v404_all_pm_un.fasta -g potato_dm_v404_all_pm_un.fasta -t 12 -s DM404


**Run ChorusHomo to design probes for close related species**

This step requires more memory than Chorus2, thus set the threads smaller will be safer.

.. code-block:: bash

    $ ChorusHomo -ga potato_dm_v404_all_pm_un.fasta -gb S_lycopersicum_chromosomes.3.00.fa \
                 -i DM404/potato_dm_v404_all_pm_un.fasta_all.bed -t 8


**Check the designed probes**

Output files will be saved to \"probes\" folder. The final file *potato_dm_v404_all_pm_un.fasta.homo.csv* 
contains the conserved probes between potato and tomato genome.

The csv file contains 10 data per line:

.. code-block:: log

    $ head probes/potato_dm_v404_all_pm_un.fasta.homo.csv
    index,probe_seq,genomeA_chr,genomeA_start,genomeA_end,genomeA_identity,genomeB_chr,genomeB_start,genomeB_end,genomeB_identity
    1,ATTTTCCATGGACCTCATTAAGATTAGCTATTGAACCAGTTACCC,chr01,97858,97902,0.99,SL3.0ch07,2272165,2272211,0.70
    2,ACAGCCAAATCGTCCCATATTCAAGGATAAACGACCCACGAATCA,chr01,103008,103052,0.99,SL3.0ch07,59525112,59525156,0.80
    3,ATCTATATCTACTACACCAGAATATTCATACACAAATAAATTACT,chr01,127095,127139,0.99,SL3.0ch06,8431503,8431548,0.83
    4,ACTGGAAGTAATATTTTGTCACAGGAATAATTCGGTAAATCTGAG,chr01,127815,127859,0.99,SL3.0ch01,27866488,27866533,0.85
    5,AGGATGAAACCGTATTCCAACCTTGCTTTTCCTTTATTGATTAAT,chr01,133524,133568,0.99,SL3.0ch10,9302271,9302317,0.79
    6,AAGCAGATATATCGTTCATCATACTTTATTTACATGGGGAAACAA,chr01,133859,133903,0.99,SL3.0ch10,39942300,39942345,0.78
    7,TTTCTTGTGCATATTTCTAAATTGTACTGTGCAAAACTTTTCCCT,chr01,134020,134064,0.99,SL3.0ch00,85946,85991,0.83
    8,AAATGATTTGCTCTTGACTGTACGTATGCCTGCCGTCTTCGTTGA,chr01,134136,134180,0.99,SL3.0ch01,63698401,63698446,0.74
    9,CCTCAAGCTTACCTACAATTAGCATAGGCAGAGTTACAAGTGGAA,chr01,134190,134234,0.99,SL3.0ch01,27863901,27863944,0.91

Column 1 is index, column 2 is probe sequence, column 3-5 are probe chromosome/start/end site in genomeA (source), 
column 6 is the identity of probe in genomeA, column 7-9 are probe chromosome/start/end site in genomeB (target), 
column 10 is the identity of probe in genomeB.

Users can select high conserved probes (for example, genomeB_identity >= 0.90) for FISH experiments in close related species.



Use ChorusNoRef to design oligo probes without a reference genome
-----------------------------------------------------------------

In this tutorial, we will design oligo probes for two wild potato species, *S. etuberosum* and *S. jamesii*, 
the two species do not have reference genomes.

Run ChorusNoRef
******************************

**Make a project folder**

.. code-block:: bash

    $ cd ~
    $ mkdir sampleproject
    $ cd sampleproject


**Download Genome file of close related species**

.. code-block:: bash

    $ wget http://solanaceae.plantbiology.msu.edu/data/potato_dm_v404_all_pm_un.fasta.zip
    $ unzip potato_dm_v404_all_pm_un.fasta.zip


**Download shotgun sequences of all species (at least 5x reads)**

.. code-block:: bash

    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/006/SRR5349606/SRR5349606_1.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/006/SRR5349606/SRR5349606_2.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/003/SRR5349573/SRR5349573_1.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/003/SRR5349573/SRR5349573_2.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/004/SRR5349574/SRR5349574_1.fastq.gz
    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR534/004/SRR5349574/SRR5349574_2.fastq.gz

SRR5349606 is from *S. tuberosum* (DM404), SRR5349573 is from *S. etuberosum*, SRR5349574 is from *S. jamesii*.


**Run Chorus2, ChorusNGSfilter and ChorusNGSselect to design probes in related species**

Run Chorus2

.. code-block:: bash

    $ Chorus2 -i potato_dm_v404_all_pm_un.fasta -g potato_dm_v404_all_pm_un.fasta -t 12

Run ChorusNGSfilter

.. code-block:: bash

    $ ChorusNGSfilter -g potato_dm_v404_all_pm_un.fasta -i SRR5349606_1.fastq.gz,SRR5349606_2.fastq.gz -t 12 \
                      -p probes/potato_dm_v404_all_pm_un.fasta_all.bed -o potato_dm_v404_all_pm_un.fasta_kmer.bed

Run ChorusNGSselect

.. code-block:: bash

    $ ChorusNGSselect -i potato_dm_v404_all_pm_un.fasta_kmer.bed -o potato_dm_v404_all_pm_un.fasta_kmerfiltered.bed


**Run ChorusNoRef to design probes in target species**

.. code-block:: bash

    $ ChorusNoRef -g potato_dm_v404_all_pm_un.fasta -p potato_dm_v404_all_pm_un.fasta_kmerfiltered.bed \
                  -r1 SRR5349573_1.fastq.gz,SRR5349574_1.fastq.gz -r2 SRR5349573_2.fastq.gz,SRR5349574_2.fastq.gz \
                  -n etuberosum,jamesii -t 12

**Check the designed probes**

Output files will be saved to \"noRefprobes\" folder. 5 files generated.

.. code-block:: log

    $ ls -lh noRefprobes
    -rw-rw-r-- 1 liu liu 5.2M 6月  28 17:14 etuberosum_indel_probe.txt
    -rw-rw-r-- 1 liu liu 104M 6月  28 17:14 etuberosum_jamesii_cns_probe.csv
    -rw-rw-r-- 1 liu liu  62M 6月  28 17:14 etuberosum_probe.txt
    -rw-rw-r-- 1 liu liu 5.2M 6月  28 17:14 jamesii_indel_probe.txt
    -rw-rw-r-- 1 liu liu  62M 6月  28 17:14 jamesii_probe.txt
    drwxrwxr-x 2 liu liu 4.0K 6月  28 16:25 tmp

*etuberosum_probe.txt* and *jamesii_probe.txt* are probes with SNPs or identical compared to DM.
*etuberosum_indel_probe.txt* and *jamesii_indel_probe.txt* are probes with indels compared to DM.
*etuberosum_jamesii_cns_probe.csv* is consensus probes among three species after quality filter.

For probes with SNPs:

.. code-block:: log

    $ head -n 3 etuberosum_probe.txt
    chr00   130544  130588  AGATTTTGCCCATTCTCATGACGCTTTTGTGATTTCAAAACTTTG   366     +       AGATTTAGCTCATTTTCATGGCGATTTTGTGATTTCAAGACTTTG   4
    chr00   129321  129365  AATACTATTAGATGATGACTAAGAGTAATGCTAGTGTATATAAAT   262     -       CTTTATATACACTAGCATTACTCTTAGTCATCATCTAATATTATT   3
    chr00   174138  174182  TTATAGTTGTCTAGGATGGAAGGGTTCTTGATTCACTGGTGTTGA   341     -       TCAACACTAGCGAATCAAGAACCCTTCCATCCTAGACAACTATAA   2

column 1-3 is location of this probe base on reference genome, column 4 is probe from reference genome, column 5 is kmer score, 
column 6 is strand, column 7 is probe for *S. etuberosum*, column 8 is how many copy can be found in etuberosum illumina reads.

For probes with indels:

.. code-block:: log

    $ head -n 3 etuberosum_probe.txt
    chr00   298161  298205  TGATGAAGGTGAAAGTAGCATAGATCATGGGGAGTTGTTTGGATT   456     +       etuberosum      TGATGAAGGTGAAAGTAGCATAGTGCATAGATCATGGGGAGTTGTTTGGATT
    chr00   298247  298291  GAATGATGAGTCAATCTGATAATTCATAGAATCAAATTTGTATGA   281     +       etuberosum      GAGTGATGAGTCAATCCATAAAGGCACCTGATAATTCATAGAATCAAATTTGTATTA
    chr00   298193  298237  TCTTTAATTTACACCATAAAGTTTACTCACAAAATCCAAACAACT   495     -       etuberosum      AGTTGTTTGGATTTTGTGAAGAGAGCAGTAAACTTTATGGTGTAAATAAAAGA

column 1-3 is location of this probe based on reference genome, column 4 is probe from reference genome, column 5 is kmer score, 
column 6 is strand, column 7 is sample name, column 8 is probe for *S. etuberosum*.

For consensus probes:

.. code-block:: log

    $ head -n 3 etuberosum_jamesii_cns_probe.csv
    chrom,start,end,refseq,etuberosum,jamesii,consensusprobe,consensusscore,consensussite,consensusdiff
    chr00,130544,130588,AGATTTTGCCCATTCTCATGACGCTTTTGTGATTTCAAAACTTTG,AGATTTAGCTCATTTTCATGGCGATTTTGTGATTTCAAGACTTTG,AGATTTAACCCATTTTCATGGCGCTTTTGTAATTTCAAGACTTTG,
    AGATTTAGCCCATTTTCATGGCGCTTTTGTGATTTCAAGACTTTG,0.9407407407407408,37,8
    chr00,129321,129365,AATACTATTAGATGATGACTAAGAGTAATGCTAGTGTATATAAAT,CTTTATATACACTAGCATTACTCTTAGTCATCATCTAATATTATT,CTTTATATACACTAGCATTACTCTTAGTCATCATCTAATATTGCT,
    CTTTATATACACTAGCATTACTCTTAGTCATCATCTAATATTAAT,0.7407407407407407,11,35

column 1-3 is location of probe based on reference genome, column 4-6 are probes in DM, etuberosum and jamesii, respectively.
consensusprobe means the consensus probe among three species.
consensusscore is calculated with fomula\: 
    (probe length * number of species - number of difference) / (probe length * number of species)
Consensussite means all identical nt in cns probe.
consensusdiff means how many nt different compare with cns probe.



Video Tutorials
-----------------

Youtube Playlist
*****************

    Playlist: https://www.youtube.com/playlist?list=PLo8q8tqFX5J27OsuKYFpd-gOtl8Qgf70X


Bilibili Playlist
*****************

    Playlist: https://www.bilibili.com/video/BV1W54y1S7qS/

