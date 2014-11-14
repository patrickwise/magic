# magic version

VERSION = 0.1

CTAGS_FILE=`gcc -M *.c -Iinclude | sed "s/ /\n/g" | grep -v "[:\\]"`

# paths
PREFIX = /home/patrick/local
MANPREFIX = ${PREFIX}/share/man
DYNPREFIX = ${PREFIX}/dynamic

# compile settings
INCS = -I. -Iinclude
CFLAGS = -std=c99 -Wall -Wextra -g ${INCS}  -pipe -march=native
LDFLAGS = -lm -g -lnlopt -latlas -lclapack
CC = gfortran

FFLAGS = -c -pipe -march=native -O2
F90 = gfortran

FSRC=$(wildcard *.f90)

SRC = ${wildcard *.c}
OBJ = ${FSRC:.f90=.o} ${SRC:.c=.o}

export CC
export CFLAGS
export LDFLAGS
export DEF

all: fortran_code options magic

options:
	@echo magic build options
	@echo "CFLAGS  = ${CFLAGS}"
	@echo "LDFLAGS = ${LDFLAGS}"
	@echo "CC      = ${CC}"
	@echo "DEF     = ${DEF}"
	@echo "OBJ     = ${OBJ}"

setup: cubature finitediff
	@echo setting up build enviroment

clean:
	rm -f ${OBJ}

dist: clean
	rm -f $(wildcard *.f90)
	rm -f $(wildcard *.mod)
	rm -f pcubature.c

fortran_code:
	$(F90) $(FSRC) $(FFLAGS)

cubature:
	@echo ">>> Downloading cubature library."
	wget "http://ab-initio.mit.edu/cubature/cubature-1.0.tgz"
	@echo ">>> Extracing library."
	tar xvzf cubature-1.0.tgz
	@echo ">>> Linking source."
	ln -sf cubature-1.0/pcubature.c
	ln -sf ../cubature-1.0/clencurt.h include
	ln -sf ../cubature-1.0/cubature.h include
	ln -sf ../cubature-1.0/vwrapper.h include
	ln -sf ../cubature-1.0/converged.h include

finitediff:
	@echo ">>> Downloading finitediff program source."
	git clone "https://github.com/bjodah/finitediff.git"
	@echo ">>> Linking source."
	ln -sf ./finitediff/finitediff/c_fornberg.f90
	ln -sf ./finitediff/finitediff/fornberg.f90
	ln -sf ../finitediff/finitediff/c_fornberg.h include

tags:
	@echo ">>> Restoring ctags."
	@ctags --fields=+l $(CTAGS_FILE)
	@sed -i "s/language:C++/language:C/g" tags #ctags defaults to c++ for header files
	@sed -i "/^__/d" tags

${OBJ}: config.mk

magic: ${OBJ}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ} ${LDFLAGS}

install:
	@mkdir -p $(PREFIX)/bin
	@cp -f magic $(PREFIX)/bin
	@mkdir -p $(DYNPREFIX)

uninstall:
	@echo uninstalling
	@rm -fr $(PREFIX)/bin/magic $(DYNPREFIX)

push:
	git add *
	git commit -m "commit"
	git push -u origin master


.PHONY: all options setup clean dist install uninstall optimize tags push dist
