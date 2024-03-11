CXX := clang++-17
CXXFLAGS := -Wall -Wextra -Werror -O2 -std=c++17
PROGNAME := main_2_1.exe
SRC := random.cpp monte_carlo.cpp utilities.cpp
OBJ := random.o monte_carlo.o utilities.o
MAIN := main_2_1.cpp 

$(PROGNAME): $(MAIN) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $^

%.o: %.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<

.PHONY: clean

clean:
	rm -r *.o *.exe