# $@ is the file name of the target
# $^ is the name of all dependencies separated by spaces, duplicate names removed
# -c means compile only (no executable, only make the .o file)a
#

CXX=g++ -std=c++11
INCLUDES = -I/users/ryansun/Documents/Research/Paper2/Software/includes
CXXFLAGS = -g -Wall $(INCLUDES)
LDFLAGS = -g -L/usr/local/lib
LDLIBS = -larmadillo -lMvtnorm

all: GOF_exact_pvalue

# 
# It doesn't automatically know the name of the binary
#

GOF_exact_pvalue: GOF_exact_pvalue.o 
	$(CXX) $(LDFLAGS) GOF_exact_pvalue.o $(LDLIBS) -o GOF_exact_pvalue 

#
# It automatically knows to make the .o from the .c, that's why we only need $^
#

GOF_exact_pvalue.o: GOF_exact_pvalue.cpp 
	$(CXX) $(CXXFLAGS) -c $^ 

clean:
	rm -f *.o GOF_exact_pvalue

