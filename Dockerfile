FROM ubuntu:latest as builder

# Installations based on saragiuliani/mum-phinder
RUN apt-get update -qq && \
    apt-get install -y zlib1g-dev \
                    git \
                    cmake \
                    build-essential \
                    python3 \
                    gcc-9 \
                    g++-9 \
                    && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 9 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9

RUN git clone https://github.com/oma219/spumoni.git; cd spumoni; mkdir build; cd build; cmake ..; make; make install;
ENV PATH=$HOME/spumoni/build:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
ENV SPUMONI_BUILD_DIR=/spumoni/build