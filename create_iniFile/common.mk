OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -g

#FFTW3DIR :=/Users/users/hutter/Libraries/include
#FFTW_CFLAGS := -I$(FFTW3DIR)
#FFTW3LIBDIR :=/Users/users/hutter/Libraries/lib
#FFTW3_LINK := -L$(FFTW3LIBDIR) -lfftw3

LDFLAGS := -lm
CFLAGS := -c -std=c99 -march=native $(WARNING) $(OPTIMIZE)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    COMPILER := clang
else
    COMPILER := gcc
endif

ifdef USE-MPI
    CC := mpicc
    CFLAGS += -D __MPI
    LDFLAGS += -lmpi
else
    CC := $(COMPILER)
    CFLAGS += 
    LDFLAGS += 
endif
