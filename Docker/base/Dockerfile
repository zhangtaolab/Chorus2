FROM forrestzhang/ubuntu-dev

MAINTAINER Tao Zhang "forrestzhang1982@gmail.com"


RUN apt-get update && apt-get install cython3


RUN mkdir /opt/software
WORKDIR /opt/software
ADD https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz /opt/software/
RUN tar zxvf jellyfish-2.2.3.tar.gz && mv jellyfish-2.2.3  jellyfish && cd jellyfish && ./configure && make && make install

#WORKDIR /opt/software
ADD https://github.com/forrestzhang/primer3-py/archive/unicode.zip /opt/software/
RUN unzip unicode.zip && cd primer3-py-unicode && python3 setup.py install

#WORKDIR /opt/software
RUN git clone https://github.com/lh3/bwa.git && cd bwa && make

RUN pip3 install numpy pyfasta

#RUN mkdir Chorus
# RUN mkdir -p Chorus/Choruslib
# ADD Chorus.py Chorus/
# ADD Choruslib/* Chorus/Choruslib/
WORKDIR /opt/software
RUN git clone https://github.com/zhangtaolab/Chorus2.git


ENV PATH /opt/software/bwa:/opt/software/jellyfish/bin:$PATH

ENV CHORUS_USER=chorus \
    CHORUS_UID=1555 \
    CHORUS_GID=1555 \
    CHORUS_HOME=/home/chorus

RUN groupadd -r $CHORUS_USER -g $CHORUS_GID
RUN useradd -u $CHORUS_UID -r -g $CHORUS_GID -d $CHORUS_HOME -c "Chorus user" $CHORUS_USER

ADD ./script/Chorus.sh /usr/bin/Chorus
RUN chmod +x /usr/bin/Chorus

WORKDIR $CHORUS_HOME


