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
WORKDIR /root
RUN git clone https://github.com/karolisr/semblans.git
WORKDIR /root/semblans
RUN ./install.sh
RUN mv bin/* /usr/local/bin/ && rm -rf bin
RUN mv lib/* /usr/local/lib/ && rm -rf lib
RUN mv external /usr/local/ && rm -rf external
RUN mv include/* /usr/local/include/ && rm -rf include
RUN cd examples && tar xzf ChloroSubSet.tar.gz && rm -f ChloroSubSet.tar.gz
RUN mkdir -p /usr/local/share/semblans-examples && mv examples/* /usr/local/share/semblans-examples && rm -rf examples
WORKDIR /usr/local/external/Rcorrector/jellyfish
RUN make clean && ./configure && make && make install
WORKDIR /root/semblans
RUN make clean && make
RUN mv bin/* /usr/local/bin/ && rm -rf bin && make clean
RUN mkdir -p /usr/local/share/semblans-ensembl-refs
RUN cp utils/semblans-download-ensembl-refs.sh /usr/local/share/semblans-ensembl-refs
WORKDIR /usr/local/share/semblans-ensembl-refs
RUN ./semblans-download-ensembl-refs.sh
WORKDIR /
ENTRYPOINT [ "semblans" ]
