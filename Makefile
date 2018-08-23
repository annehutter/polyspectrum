SOURCES := 	./src/main.c \
		./src/domain.c \
		./src/xmem.c \
		./src/xstring.c \
		./src/parse_ini.c \
		./src/utils.c\
		./src/confObj.c \
		./src/grid.c \
		./src/fft.c \
		./src/filter.c \
		./src/utils_fftw.c \
		./src/polyspectrum.c


OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := polyspectrum

USE-MPI=YES
 
include common.mk


.PHONY: all clean clena celan celna

all: $(SOURCES) $(EXECUTABLE)

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
