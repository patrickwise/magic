# magic version

VERSION = 0.1

CTAGS_FILE=`gcc -M *.c -Iinclude | sed "s/ /\n/g" | grep -v "[:\\]"`

# paths
PREFIX = /home/patrick/local
MANPREFIX = ${PREFIX}/share/man
DYNPREFIX = ${PREFIX}/dynamic

# compile settings
INCS = -I. -Iinclude
CFLAGS = -std=c99 -Wall -Wextra -g -pedantic ${INCS}  -pipe -march=native
LDFLAGS = -lm -lnlopt -latlas -g
CC = gfortran

FFLAGS = -c -pipe -march=native -O2
F90 = gfortran

FSRC=$(wildcard *.f90)

SRC = ${wildcard *.c}
OBJ = ${FSRC:.f90=.o} ${SRC:.c=.o}

DEF = -DDYNPREFIX=\"$(DYNPREFIX)/\"

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
#build .mod files. This will fail the first time through.
	@gfortran ${FSRC} 2> /dev/null || true

clean:
	rm -f ${OBJ}

dist: clean
	rm -f $(wildcard *.f90)
	rm -f $(wildcard *.mod)
	rm -f $(wildcard *.tgz)
	rm -f pcubature.c
	rm -rf finitediff cubature-1.0
	rm -rf $(wildcard dynamic/*)

fortran_code:
	@$(F90) $(FSRC) $(FFLAGS) || true #|| echo "Did you run 'make setup'?" && false

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

.c.o:
	${CC} ${CFLAGS} ${DEF} ${INC} $< -c -o $@

magic: ${OBJ}
	@echo CC -o $@
	@${CC} -o $@ ${OBJ} ${LDFLAGS} ${DEF}

install: all
	@mkdir -p $(PREFIX)/bin
	@cp -f magic plotter.sh $(PREFIX)/bin
	@chmod 755 ${PREFIX}/bin/{magic,plotter.sh}
	@mkdir -p $(DYNPREFIX)
	@chmod 755 ${DYNPREFIX}

uninstall:
	@echo uninstalling
	@rm -fr $(PREFIX)/bin/magic $(DYNPREFIX)

push: dist
	git add *
	git commit
	git push -u origin master


.PHONY: all options setup clean dist install uninstall optimize tags push dist
