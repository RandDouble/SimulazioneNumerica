PROGNAME= main_1_2.exe
MAIN= main_1_2.cpp 
SRC_FILES = random.cpp utilities.cpp

# Compiler and compiler flags
# Default compiler is Clang++-17 with libc++ as standard library implementation
INCFLAGS= -I../../src
CXXFLAGS= $(INCFLAGS) -Wall -Wextra -Werror -std=c++17 -g3 -O2 -march=native -ffast-math -DUSE_PARALLEL
# O2 because it can show more warning due to flow variable analisys, without taking too much compilation time.
# c++17 because it has some niceties, like in better algorithm and numeric, enabling parallelism.

ifeq ($(COMPILER), CLANG_LIBCPP) # Only because it support c++20 needed for other projects
	CXX = clang++-17
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CXX = clang++-17
	CXXFLAGS += -stdlib=libstdc++
	LDFLAGS += -ltbb
else ifeq ($(COMPILER), GNU)
	CXX = g++
	LDFLAGS += -ltbb
else # Defaults to CLANG_LIBCPP, flags to set library and use parallel algorithms
	CXX = clang++-17
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
endif


SOURCE_DIR= ../../src
OBJ_DIR= ./obj
SRC= $(addprefix $(SOURCE_DIR), $(SRC_FILES))
OBJ= $(addprefix $(OBJ_DIR)/,$(addsuffix .o , $(basename $(SRC_FILES))))

.PHONY: all
all : $(PROGNAME)

$(PROGNAME): $(MAIN) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) -c -o $@ $(CXXFLAGS) $<

.PHONY: clean

clean:
	rm -r $(OBJ) *.exe