FROM nvidia/cuda:11.0.3-devel-centos8 as BUILDER

# Preparing the environment
RUN yum clean all && yum --enablerepo=extras install -y epel-release && yum update -y
RUN yum install -y git wget zlib-devel boost boost-devel cmake3 openssl openssl-devel
RUN yum groupinstall "Development Tools" -y   
RUN yum clean all
RUN rm -rf /var/cache/yum/*

# Installing LASTZ (SegAlign's dependency)
WORKDIR /usr/local/src/
RUN wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
RUN tar -zxf lastz-1.04.03.tar.gz && rm lastz-1.04.03.tar.gz
RUN make -C lastz-distrib-1.04.03/src
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
## Cloning it
WORKDIR /usr/local/src/
RUN git clone https://github.com/thiagogenez/SegAlign.git --recursive 
WORKDIR SegAlign

## Download TBB (Segalign's build dependency)
RUN wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_lin.tgz
RUN tar -xvf tbb2019_20191006oss_lin.tgz && rm -f tbb2019_20191006oss_lin.tgz

## Building Segalign
WORKDIR build
RUN cmake3 -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake ..
RUN make 
RUN cp segalign* /usr/local/bin 
RUN cp ../scripts/run_segalign* /usr/local/bin
RUN cp /usr/local/src/SegAlign/tbb*/lib/intel64/gcc*/lib* /usr/local/lib/

# Create a thinner final Docker image with only runtime dependencies
# FROM nvidia/cuda:10.2-runtime-centos7


# Install runtime dependencies
#RUN yum clean all && yum update -y
#RUN yum install -y make zlib wget boost boost-devel krb5-libs
#RUN yum groupinstall "Development Tools" -y   
#RUN yum clean all
#RUN rm -rf /var/cache/yum/*

# Installing Paralell (SegAlign's runtime dependency)
WORKDIR /usr/local/src/
RUN wget https://ftp.gnu.org/gnu/parallel/parallel-20201122.tar.bz2
RUN tar xjf parallel-20201122.tar.bz2 && rm parallel-20201122.tar.bz2
WORKDIR parallel-20201122
RUN ./configure && make && make install

# copy all the compiled binaries from BUILDER
#COPY --from=BUILDER /usr/local/bin /usr/local/bin

# copy the tbb shared library
#COPY --from=BUILDER /usr/local/src/SegAlign/tbb*/lib/intel64/gcc*/lib* /usr/local/lib/

# add the library path
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

# Go to workspace
WORKDIR /workspace
