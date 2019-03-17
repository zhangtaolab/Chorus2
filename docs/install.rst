Install
=======

Docker
------

Install Docker_

.. _Docker: https://docs.docker.com/engine/installation/

Terminal version
*****************

Download Chorus2:

.. code-block:: bash

    $ docker pull forrestzhang/docker-chorus

    $ docker run -v $PWD:/home/chorus -e CHORUS_USER=$USER -e CHORUS_UID=$UID forrestzhang/docker-chorus -h

    usage: Chorus [-h] [--version] [-j JELLYFISH] [-b BWA] -g GENOME -i INPUT
              [-s SAVED] [-p PRIMER] [-t THREADS] [-l LENGTH]
              [--homology HOMOLOGY] [-d DTM] [--docker DOCKER]

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
    --docker DOCKER

GUI (Test Version)
*******************

Download Chorus2:

.. code-block:: bash

    $ docker pull forrestzhang/chorus-gui

    $ docker run -i -t -p 6080:80 -v $PWD:/root/Desktop/Data forrestzhang/chorus-gui

Open your web browse, enter http://localhost:6080

In your web browse, open a LXTerminal

.. code-block:: bash

    $ python3 /opt/software/Chorus/ChorusGUI.py

.. image:: _static/docker_GUI.jpg

Manually Install
----------------

Ubuntu 14.04 (terminal)
***********************

Install dependent package

.. code-block:: bash

    $ apt-get update && apt-get install  -y build-essential \
        cython3 \
        zlib1g-dev \
        zlibc \
        openjdk-7-jre \
        git \
        libboost-dev \
        autoconf \
        libncursesw5-dev \
        libncurses5 \
        ncurses-dev \
        libboost-thread-dev \
        python3-pip \
        samtools \
        unzip \
        python \
        curl \
        wget


Install jellyfish

.. code-block:: bash

    $ mkdir /opt/software

    $ cd /opt/software

    $ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz

    $ tar zxvf jellyfish-2.2.3.tar.gz

    $ mv jellyfish-2.2.3  jellyfish

    $ cd jellyfish

    $ ./configure && make && make install


Install bwa

.. code-block:: bash

    $ cd /opt/software

    $ git clone https://github.com/lh3/bwa.git

    $ cd bwa

    $ make


Install primer3-py

.. code-block:: bash

    $ cd /opt/software

    $ wget https://github.com/forrestzhang/primer3-py/archive/unicode.zip

    $ unzip unicode.zip

    $ cd primer3-py-unicode

    $ python3 setup.py install

Download Chorus2

.. code-block:: bash

    $ cd /opt/software

    $ git clone https://github.com/zhangtaolab/Chorus2.git

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

Anaconda (Terminal)
********************************

Download and Install Anaconda_ (python 3.x verison)

.. _Anaconda: https://www.continuum.io/downloads


Add bioconda

.. code-block:: bash

    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda


Install dependent package

.. code-block:: bash

    conda install bwa
    conda install jellyfish
    pip install Cython 
    pip install primer3-py
    pip install pyfasta
    pip install pyBigWig


Ubuntu 14.04 (GUI) Test Version
********************************

Install dependent package

.. code-block:: bash

    $ apt-get update && apt-get install -y cython3  build-essential \
        zlib1g-dev \
        zlibc \
        git \
        libboost-dev \
        autoconf \
        libncursesw5-dev \
        libncurses5 \
        ncurses-dev \
        libboost-thread-dev \
        python3-pip \
        samtools \
        unzip \
        python \
        curl \
        wget \
        python3-pyqt5 \
        libfreetype6-dev \
        libxft-dev \
        python3-matplotlib

    $ apt-get remove -y python3-matplotlib


Install jellyfish

.. code-block:: bash

    $ mkdir /opt/software

    $ cd /opt/software

    $ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz

    $ tar zxvf jellyfish-2.2.3.tar.gz

    $ mv jellyfish-2.2.3  jellyfish

    $ cd jellyfish

    $ ./configure && make && make install


Install bwa

.. code-block:: bash

    $ cd /opt/software

    $ git clone https://github.com/lh3/bwa.git

    $ cd bwa

    $ make


Install primer3-py

.. code-block:: bash

    $ cd /opt/software

    $ wget https://github.com/forrestzhang/primer3-py/archive/unicode.zip

    $ unzip unicode.zip

    $ cd primer3-py-unicode

    $ python3 setup.py install


Install Python dependent package

.. code-block:: bash

    $ pip3 install numpy pyfasta matplotlib

    $ pip3 install pandas

Download and run ChorusGUI

.. code-block:: bash

        $ cd /opt/software

        $ git clone https://github.com/zhangtaolab/Chorus2.git

        $ python3 /opt/software/Chorus2/ChorusGUI.py