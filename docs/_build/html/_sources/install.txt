Install
=======

Docker
------

Install Docker_

.. _Docker: https://docs.docker.com/engine/installation/

Terminal version
*****************

Download Chorus:

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

Download Chorus:

.. code-block:: bash

    $ docker pull forrestzhang/chorus-gui

    $ docker run -i -t -p 6080:6080 -v $PWD:/home/ubuntu/Data forrestzhang/chorus-gui

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

Download Chorus

.. code-block:: bash

    $ cd /opt/software

    $ git clone https://github.com/forrestzhang/Chorus.git

    $ python3 python3 /opt/software/Chorus/Chorus.py -h



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

    $ pip3 install pandas==0.16.2

Download and run Chorus

.. code-block:: bash

        $ cd /opt/software

        $ git clone https://github.com/forrestzhang/Chorus.git

        $ python3 /opt/software/Chorus/ChorusGUI.py