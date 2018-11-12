CXX=g++-4.9
CXXFLAGS=-std=c++11 -O3

MIPS_OBJ=$(patsubst %.cc,%.o,$(wildcard src/mips/*.cc))
BRANCHBOUND_OBJ=$(patsubst %.cc,%.o,$(wildcard src/branchbound/*.cc))
zetaIL_DIR=reproduce_results/src/zetaIL
SPP_DIR=reproduce_results/src/SPP

all: src/train_WHInter $(zetaIL_DIR)/train_zetaIL $(SPP_DIR)/train_SPPbreadthFirst

src/train_WHInter: src/train_WHInter.o $(MIPS_OBJ) $(BRANCHBOUND_OBJ) src/other/func.o
	$(CXX) -o src/train_WHInter src/train_WHInter.o src/other/func.o src/branchbound/*.o src/mips/*.o

src/train_WHInter.o: src/train_WHInter.cc
	$(CXX) -c $(CXXFLAGS) -o src/train_WHInter.o src/train_WHInter.cc

src/mips/%.o: src/mips/%.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

src/branchbound/%.o: src/branchbound/%.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

src/other/func.o: src/other/func.cc
	$(CXX) -c $(CXXFLAGS) -o src/other/func.o src/other/func.cc


$(zetaIL_DIR)/train_zetaIL: $(zetaIL_DIR)/train_zetaIL.o $(MIPS_OBJ) src/other/func.o
	$(CXX) -o $(zetaIL_DIR)/train_zetaIL $(zetaIL_DIR)/train_zetaIL.o src/other/func.o src/mips/*.o
$(zetaIL_DIR)/train_zetaIL.o: $(zetaIL_DIR)/train_zetaIL.cc
	$(CXX) -c $(CXXFLAGS) -o $(zetaIL_DIR)/train_zetaIL.o $(zetaIL_DIR)/train_zetaIL.cc


$(SPP_DIR)/train_SPPbreadthFirst: $(SPP_DIR)/train_SPPbreadthFirst.cc
	$(CXX) $(CXXFLAGS) -o $(SPP_DIR)/train_SPPbreadthFirst $(SPP_DIR)/train_SPPbreadthFirst.cc

clean:
	find . -type f -name '*.o' -exec rm {} \;
	#-rm -f src/*.o
