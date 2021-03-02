FROM nvidia/cuda:11.2.1-devel-centos7 as  builder

# Preparing the environment
RUN yum update 
RUN yum install -y make sudo git wget curl file perl-core zlib-devel boost boost-devel 
RUN yum groupinstall "Development Tools" -y   
RUN yum clean all
RUN rm -rf /var/cache/yum/*

# Installing OpenSSL (cmake's dependency)
WORKDIR /usr/local/src/
RUN wget https://www.openssl.org/source/openssl-1.1.1j.tar.gz
RUN tar -xf openssl-1.1.1j.tar.gz && rm openssl-1.1.1j.tar.gz
WORKDIR openssl-1.1.1j
RUN ./config --prefix=/usr --openssldir=/etc/ssl --libdir=lib no-shared zlib-dynamic
RUN make && make test && make install

# Installing cmake (SegAlign's dependency)
WORKDIR /usr/local/src/
RUN wget https://github.com/Kitware/CMake/releases/download/v3.19.6/cmake-3.19.6.tar.gz
RUN tar zxvf cmake-3.19.6.tar.gz && rm cmake-3.19.6.tar.gz
WORKDIR cmake-3.19.6
RUN ./bootstrap --prefix=/usr/local && sudo make && sudo make install

# Installing LASTZ (SegAlign's dependency)
WORKDIR /usr/local/src/
RUN wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
RUN tar -zxf lastz-1.04.03.tar.gz && rm lastz-1.04.03.tar.gz
RUN make -C lastz-distrib-1.04.03/src
RUN ls lastz-distrib-1.04.03/src
RUN cp lastz-distrib-1.04.03/src/lastz /usr/local/bin
RUN rm -rf lastz-distrib-1.04.03

# Installing kentUtils (faToTwoBit and twoBitToFa) (SegAlign's dependency) 
WORKDIR /usr/local/src/
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faToTwoBit
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa
RUN chmod +x faToTwoBit
RUN chmod +x twoBitToFa
RUN mv faToTwoBit /usr/local/bin/
RUN mv twoBitToFa /usr/local/bin/


# Installing Segalign
## Cloning
WORKDIR /usr/local/src/
RUN git clone https://github.com/thiagogenez/SegAlign.git --recursive 

## Change the branch to the one pointed by the commit
WORKDIR SegAlign
COPY segalign.commit /
RUN git checkout $(cat /segalign.commit)
RUN git rev-parse HEAD > /segalign.commit

## Download TBB
RUN wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_lin.tgz
RUN tar -xvf tbb2019_20191006oss_lin.tgz && rm tbb2019_20191006oss_lin.tgz

## Building Segalign
WORKDIR build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake ..
RUN make 
RUN cp segalign* /usr/local/bin 
RUN cp ../scripts/run_segalign* /usr/local/bin



# Create a thinner final Docker image with only runtime dependencies
#FROM nvidia/cuda:11.2.1-runtime-ubuntu20.04

#RUN DEBIAN_FRONTEND=noninteractive apt install -y tzdata

# Install runtime dependencies
#RUN DEBIAN_FRONTEND=noninteractive \
#    TZ=Europe/London \
#    apt-get -qq -y update && \
#    apt-get -qq -y upgrade && \
#    apt-get -qq -y install \
#    libkrb5-3 \
#    libk5crypto3 \
#    libboost-dev \
#    libboost-program-options-dev \
#    zlib1g \
#    parallel

# copy all the binaries
#COPY --from=builder /usr/local/bin /usr/local/bin

# copy the tbb shared library
#COPY --from=builder /WGA_GPU/tbb*/lib/intel64/gcc*/lib* /usr/local/lib/

# add the library path
#ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

# remember that commit
#COPY --from=builder /segalign.commit /

# Go to workspace
RUN mkdir -p /workspace
WORKDIR /workspace
