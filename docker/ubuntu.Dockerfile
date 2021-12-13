FROM ubuntu:20.04

RUN apt-get update -qq \
&& DEBIAN_FRONTEND=noninteractive apt-get install -yq \
 git wget build-essential cmake mercurial \
 # Comment to try with gcc/g++-9
 gcc-10 g++-10 \
&& update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 10 \
&& update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 10 \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Create lib directory
WORKDIR /home/lib
COPY . .
RUN cmake -S . -Bbuild
RUN cmake --build build --target all
RUN cd build && ctest --verbose
RUN cmake --build build --target install -v -- DESTDIR=install
