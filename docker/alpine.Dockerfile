FROM alpine:edge

RUN apk add --no-cache git build-base cmake mercurial

# Create lib directory
WORKDIR /home/lib
COPY . .

RUN cmake -S . -Bbuild
RUN cmake --build build --target all
RUN cd build && ctest --verbose
RUN cmake --build build --target install -v -- DESTDIR=install
