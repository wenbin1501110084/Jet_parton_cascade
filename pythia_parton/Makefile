EXE=mymain06
GCC=g++

PYTHIAINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --includedir}
PYTHIALIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/Berkeley_work/v2_in_jet/parton_shower_in_pythia8/pythia8310/bin/pythia8-config --libdir}

FASTJETINC=${shell /home/wenbin/Downloads/Wenbin_working/Work/code/fastjet_install/bin/fastjet-config --cxxflags}
FASTJETLIBS=${shell /home/wenbin/Downloads/Wenbin_working/Work/code/fastjet_install/bin/fastjet-config --libs}


all: $(EXE)

mymain06.o: mymain06.cc
	$(GCC) -I$(PYTHIAINC) $(FASTJETINC) -c mymain06.cc

$(EXE): mymain06.o
	$(GCC) -O3 -I$(PYTHIAINC) $(FASTJETINC) -o $(EXE) mymain06.o $(PYTHIALIBS)/libpythia8.a -ldl $(FASTJETLIBS) 

.PHONY:clean
clean:
	-rm *.o $(EXE)
