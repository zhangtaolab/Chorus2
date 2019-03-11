Tutorial
========

In this tutorial, we will build oligo probe set for Arabidopsis genome.

Using Docker Terminal Version
-----------------------------

Install
*******

Install Docker_

.. _Docker: https://docs.docker.com/engine/installation/

Download Chorus2:

.. code-block:: bash

    $ docker pull forrestzhang/docker-chorus

Parameter of Chorus2:

.. code-block:: bash

    -g GENOME, --genome GENOME
                            fasta format genome file
    -i INPUT, --input INPUT
                            fasta format input file
    -s SAVED, --save SAVED
                            result saved folder
    -p PRIMER, --primer PRIMER
                            5' labeled R primer
    -t THREADS, --threads THREADS
                            threads number or how may cpu you wanna use
    -l LENGTH, --length LENGTH
                            probe length
    --homology HOMOLOGY   homology, from 50 to 100
    -d DTM, --dtm DTM     dTm, from 0 to 37

Download Reference Genome file:

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


Using Manually Install Version
------------------------------

Install
*******

* `Manually Install`_
.. _`Manually Install`: \install.html#ubuntu-14-04-terminal

Run In Terminal
***************

Make a project folder

.. code-block:: bash

    $ cd ~
    $ mkdir sampleproject
    $ cd sampleproject

Download reference genome

.. code-block:: bash

    $ wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas


Test chorus2 software

.. code-block:: bash

    $ python3 /opt/software/Chorus2/Chorus.py -h

        usage: Chorus [-h] [--version] [-j JELLYFISH] [-b BWA] -g GENOME -i INPUT
                      [-s SAVED] [-p PRIMER] [-t THREADS] [-l LENGTH]
                      [--homology HOMOLOGY] [-d DTM] [--step STEP] [--docker DOCKER]
                      [--ploidy PLOIDY]

        Chorus Software for Oligo FISH probe design

        optional arguments:
          -h, --help            show this help message and exit
          --version             show program's version number and exit
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
                                A specific 5' labeled R primer for PCR reaction. For
                                example: CGTGGTCGCGTCTCA. (Default is none)
          -t THREADS, --threads THREADS
                                Number of threads or CPUs to use. (Default: 1)
          -l LENGTH, --length LENGTH
                                The probe length. (Default: 45)
          --homology HOMOLOGY   The maximum homology(%) between target sequence and
                                probe, range from 50 to 100. (Default: 75)
          -d DTM, --dtm DTM     The minimum value of dTm (hybrid Tm - hairpin Tm),
                                range from 0 to 37. (Default: 10)
          --step STEP           The step length for k-mer searching in a sliding
                                window, step length>=1. (Default: 5)
          --docker DOCKER       Only used in Docker version of Chorus
          --ploidy PLOIDY       The ploidy of the given genome (test version).
                                (Default: 2)

        Example:
          python3 Chorus.py -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 4 \
                            -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa -s sample

Run chorus2 software

.. code-block:: bash

    $ python3 /opt/software/Chorus/Chorus.py -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 12 \
      -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa

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

Further filter using ChorusNGSfilter

Before running ChorusNGSfilter, a set of whole-genome shotgun sequencing data is required. Here we download the shotgun reads of Arabidopsis with the accession number SRR5658649.

.. code-block:: bash

    $ wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR565/009/SRR5658649/SRR5658649_1.fastq.gz
    $ wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR565/009/SRR5658649/SRR5658649_2.fastq.gz

    $ python3 /opt/software/Chorus/ChorusNGSfilter.py -i SRR5658649_1.fq.gz,SRR5658649_2.fq.gz -z gz \
      -g TAIR10_chr_all.fas -t 12 -j /opt/software/jellyfish/bin/jellyfish \
      -p probes/TAIR10_chr_all.fas_all.bed -o probes/TAIR10_chr_all_SRR5658649.bed

After running NGS filtering, three files (*.jf, *.bw, *.bed) will output to working directory.

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

Automatic probe selection using ChorusNGSselect

.. code-block:: bash

    $ python3 /opt/software/Chorus/ChorusNGSselect.py -i probes/TAIR10_chr_all_SRR5658649.bed \
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