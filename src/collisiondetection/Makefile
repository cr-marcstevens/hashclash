CXX ?= gcc
CXXFLAGS ?= -O3 -march=native

detectcoll: detectcoll.c detectcoll.h
	$(CXX) $(CXXFLAGS) detectcoll.c -o detectcoll

clean:
	rm detectcoll
