-include local.mk

CFLAGS += -std=c99 -Wall -Wextra -Werror -Wno-unknown-pragmas -pedantic
LDFLAGS += 
LDLIBS += -lm

ifdef DEBUG
CFLAGS += -O0 -g -DDEBUG
else
CFLAGS += -Ofast -march=native -mfpmath=sse
endif

ifdef OPENMP
CFLAGS += -fopenmp
endif

.PHONY: all clean

all: nzz

clean:
	$(RM) nzz

nzz: nzz.c io.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $< $(LDLIBS)
