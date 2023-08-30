CXX=mpic++
# CCFLAGS=-Wall -g
# activate for compiler optimizations:
CCFLAGS=-Wall --std=c++17 -O3 -g -Wall
LDFLAGS=

all: life

.PHONY : tiny_test

life: life.cpp
	$(CXX) $(CCFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o life

tiny_test: life 5_example.txt 3_example.txt
	mpirun -np 1 ./life 5 5 1 5_example.txt output.txt
	@cmp --silent output.txt 5_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS TINY TEST"; else echo "Fail tiny test"; fi

	mpirun -np 1 ./life 3 3 1 3_example.txt output.txt
	@cmp --silent output.txt 3_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS TINY TEST"; else echo "Fail tiny test"; fi



