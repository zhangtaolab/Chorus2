FROM dorowu/ubuntu-desktop-lxde-vnc


MAINTAINER Tao Zhang "forrestzhang1982@gmail.com"

RUN apt-get update && apt-get install -y  wget bzip2

RUN mkdir /opt/download

ADD https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh /opt/download

WORKDIR /opt/download

RUN bash /opt/download/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm /opt/download/Miniconda3-latest-Linux-x86_64.sh

RUN mkdir -p /home/ubuntu/Desktop

RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc && \
    echo "conda activate base" >> /root/.bashrc &&\
    /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda update -y -n base conda && \
    /opt/conda/bin/conda install -y matplotlib pyfasta pysam pyqt jellyfish bwa pandas pip primer3-py git

# RUN apt-get install -y build-essential git

# RUN  /opt/conda/bin/pip install primer3-py

ENV PATH $PATH:/opt/conda/bin
RUN mkdir /opt/software
WORKDIR /opt/software
RUN git clone https://github.com/zhangtaolab/Chorus2.git

RUN echo "/opt/conda/bin/python /opt/software/Chorus2/ChorusGUI.py" > /home/ubuntu/Desktop/ChorusGUI.sh
RUN echo "/opt/conda/bin/python /opt/software/Chorus2/ChorusPBselect.py" > /home/ubuntu/Desktop/ChorusPBselect.sh
RUN chmod +x /home/ubuntu/Desktop/ChorusGUI.sh /home/ubuntu/Desktop/ChorusPBselect.sh

RUN mkdir /root/Desktop
RUN echo "/opt/conda/bin/python /opt/software/Chorus2/ChorusGUI.py" > /root/Desktop/ChorusGUI.sh
RUN echo "/opt/conda/bin/python /opt/software/Chorus2/ChorusPBselect.py" > /root/Desktop/ChorusPBselect.sh
RUN chmod +x /root/Desktop/ChorusGUI.sh /root/Desktop/ChorusPBselect.sh


VOLUME /root/Desktop/Data
ENTRYPOINT ["/startup.sh"]
