EXE=mymain07
GCC=g++

PYTHIAINC=${shell /wsu/home/he/he92/he9215/XETSCAPE/code_running/software_code/pythia8306_installed/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /wsu/home/he/he92/he9215/XETSCAPE/code_running/software_code/pythia8306_installed/bin/pythia8-config --libdir}

FASTJETINC=${shell /wsu/home/he/he92/he9215/XETSCAPE/code_running/software_code/fastjet/fastjet_install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell /wsu/home/he/he92/he9215/XETSCAPE/code_running/software_code/fastjet/fastjet_install/bin/fastjet-config --libs}


all: $(EXE)

mymain07.o: mymain07.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -c mymain07.cc

$(EXE): mymain07.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC) -o $(EXE) mymain07.o $(PYTHIALIBS)/libpythia8.a -ldl $(FASTJETLIBS) 

.PHONY:clean
clean:
	-rm *.o $(EXE)
