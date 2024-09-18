FROM ubuntu:22.04
RUN apt-get update
RUN apt-get full-upgrade
RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get -o Dpkg::Options::="--force-confold" install --yes \
--allow-downgrades --allow-remove-essential --allow-change-held-packages \
bowtie2 \
build-essential \
cmake \
default-jre \
git \
libcurlpp-dev \
librapidxml-dev \
python3 \
python3-dev \
python3-numpy \
python3-pip \
r-base \
samtools \
unzip \
wget
RUN cd /root && git clone https://github.com/karolisr/semblans.git
RUN cd /root/semblans && ./install.sh
RUN cd /root/semblans && mv bin/* /usr/local/bin/ && rm -rf bin
RUN cd /root/semblans && mv lib/* /usr/local/lib/ && rm -rf lib
RUN cd /root/semblans && mv external /usr/local/ && rm -rf external
RUN cd /root/semblans && mv include/* /usr/local/include/ && rm -rf include
RUN cd /root/semblans/examples && tar xzf ChloroSubSet.tar.gz && rm -f ChloroSubSet.tar.gz
RUN mkdir -p /usr/local/share/semblans-examples && cd /root/semblans && mv examples/* /usr/local/share/semblans-examples && rm -rf examples
RUN cd /usr/local/external/Rcorrector/jellyfish && make clean && ./configure && make && make install
RUN cd /root/semblans && make clean && make
RUN cd /root/semblans && mv bin/* /usr/local/bin/ && rm -rf bin && make clean
ENTRYPOINT [ "semblans" ]
