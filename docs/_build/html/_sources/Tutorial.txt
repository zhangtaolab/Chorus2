Tutorial
========

In this tutorial, we will build oligo porbe set for Arabidopsis genome.

Using Docker Terminal Version
-----------------------------

Install
*******

Install Docker_

.. _Docker: https://docs.docker.com/engine/installation/

Download Chorus:

.. code-block:: bash

    $ docker pull forrestzhang/docker-chorus

Parameter of Chorus:

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

*TAIR10_chr_all.fas.bed* is the probe file.

.. code-block:: log

    $ more probes/TAIR10_chr_all.fas.bed
    1      	52     	96     	TCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTA
    1      	211    	255    	TTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTC
    1      	346    	390    	CCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACT
    1      	426    	470    	TTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAAT
    1      	496    	540    	TCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGG
    1      	551    	595    	TAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAG

There are four columns in a row, first column is chromosome name, second is oligo start site, third is oligo end site, last one is oligo probe sequence. You can use excel or text editor to open this file.


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


Test chorus software

.. code-block:: bassh

    $ python3 /opt/software/Chorus/Chorus.py -h 

        usage: Chorus [-h] [--version] [-j JELLYFISH] [-b BWA] -g GENOME -i INPUT
                [-s SAVED] [-p PRIMER] [-t THREADS] [-l LENGTH]
                [--homology HOMOLOGY] [-d DTM] [--step STEP] [--docker DOCKER]

        Chorus Software for Oligo FISH probe design

        optional arguments:
        -h, --help            show this help message and exit
        --version             show program's version number and exit
        -j JELLYFISH, --jellyfish JELLYFISH
                                jellyfish path
        -b BWA, --bwa BWA     bwa path
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
        --step STEP           step length, min=1
        --docker DOCKER

Run chorus software

.. code-block:: bash

    $ python3 /opt/software/Chorus/Chorus.py -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 12 \
      -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa -s sample

When job finish, the oligo probes will output to 'sample' folder 

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

*TAIR10_chr_all.fas.bed* is the probe file.

.. code-block:: log

    $ more probes/TAIR10_chr_all.fas.bed
    1      	52     	96     	TCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTA
    1      	211    	255    	TTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTC
    1      	346    	390    	CCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACT
    1      	426    	470    	TTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAAT
    1      	496    	540    	TCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGG
    1      	551    	595    	TAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAG

There are four columns in a row, first column is chromosome name, second is oligo start site, third is oligo end site, last one is oligo probe sequence. You can use excel or text editor to open this file.
