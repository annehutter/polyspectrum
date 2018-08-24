OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -g

FFTW3DIR :=/usr/local/include
FFTW_CFLAGS := -I$(FFTW3DIR)
FFTW3LIBDIR :=/usr/local/lib64
FFTW3_LINK := -L$(FFTW3LIBDIR) -lfftw3

LDFLAGS := -lm $(FFTW3_LINK)
CFLAGS := -c -std=c99 -march=native $(WARNING) $(OPTIMIZE) $(FFTW_CFLAGS) 

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    COMPILER := clang
else
    COMPILER := gcc
endif

ifdef USE-MPI
    CC := mpicc
    CFLAGS += -D __MPI
    LDFLAGS += -lmpi -lfftw3_mpi
else
    CC := $(COMPILER)
    CFLAGS += 
    LDFLAGS += 
endif
