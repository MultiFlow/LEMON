FROM archlinux:latest

RUN pacman -Syu --noconfirm git base-devel cmake

# Create lib directory
WORKDIR /home/lib
COPY . .

RUN cmake -S . -Bbuild
RUN cmake --build build --target all
RUN cd build && ctest --verbose
RUN cmake --build build --target install -v -- DESTDIR=install
