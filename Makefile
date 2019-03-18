
CXX=mpicxx

all: VantagePointAllK.o
	$(CXX) VantagePointAllK.cpp -o vp

clean:
	rm -f *.o *~ vp
