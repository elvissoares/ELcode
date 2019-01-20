# Specify the source files, the target files, 
# and the install directory 
DIR = ELcode
SOURCES = delayed.cpp star.cpp eos.cpp eos/fermions.cpp eos/fermi.cpp thermonuclear.cpp conductivity.cpp
OUTPUTFILE  = delayed.exe
INSTALLDIR  = bin

#CPPUNIT_PATH = opt/local
CPPUNIT_PATH = usr/local
INCLUDE_DIR = -I/$(CPPUNIT_PATH)/include
LIB_DIR = -L/$(CPPUNIT_PATH)/lib
LIBS = -fopenmp

CC = g++
LFLAGS = -o

# This is the default target, which will be built when 
# you invoke make
.PHONY: all
all: $(OUTPUTFILE)

# $(SOURCES:.cpp=.o)
# This rule tells make how to build 'OUTPUTFILE' from 'SOURCES'
$(OUTPUTFILE):  $(SOURCES)
	$(CC) $(LFLAGS) $@ $^ $(LIBS)

# This rule tells make to copy OUTPUTFILE to the binaries subdirectory,
# creating it if necessary
.PHONY: install
install:
	mkdir -p $(INSTALLDIR)
	mv $(OUTPUTFILE) $(INSTALLDIR)

# This rule tells make to run OUTPUTFILE in the binaries subdirectory,
.PHONY: clean
clean:
	rm $(OUTPUTFILE)
