# @author: Stefano Pilosio
# @scope: Laboratorio di Simulazione Numerica 2023 / 2024

PROGNAME= main_9_1.exe
MAIN= main_9_1.cpp
SRC_FILES = random.cpp utilities.cpp
HEADER_FILES =
LOCAL_SRC_FILE =
LOCAL_HEADER_FILE = initializer.h mutation.h population.h swap_functions.h

# Compiler and compiler flags
# Default compiler is Clang++-17 with libc++ as standard library implementation
INCFLAGS= -I../../src
CXXFLAGS= $(INCFLAGS) -Wall -Wextra -Werror -O3 -ffast-math -march=native -std=c++17 -fopenmp -DUSE_PARALLEL
LDFLAGS = -larmadillo
# O2 because it can show more warning due to flow variable analisys
# c++17 because it has some niceties, like in better algorithm and numeric, enabling parallelism

ifeq ($(COMPILER), CLANG_LIBCPP) # Only because it support c++20 needed for other projects
	CXX = clang++-17
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CXX = clang++-17
	CXXFLAGS += -stdlib=libstc++
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

SRC= $(addprefix $(SOURCE_DIR)/, $(SRC_FILES))
HEAD= $(addprefix $(SOURCE_DIR)/, $(HEADER_FILES))
HEAD += $(LOCAL_HEADER_FILE)

OBJ= $(addprefix $(OBJ_DIR)/,$(addsuffix .o , $(basename $(SRC_FILES))))
SRC += $(LOCAL_SRC_FILE)
OBJ += $(addprefix $(OBJ_DIR)/, $(addsuffix .o , $(basename $(LOCAL_SRC_FILE))))


.PHONY: all
all : $(PROGNAME)

test: test_permutations.exe

test_permutations.exe: test_permutations.cpp $(HEAD) $(OBJ)
	$(CXX) -o $@ $< -g -Wall -Wextra -std=c++17 -stdlib=libc++ -fexperimental-library -I../../src $(LDFLAGS) $(OBJ)

$(PROGNAME): $(MAIN) $(OBJ) $(HEAD)
	$(CXX) -o $@ $(CXXFLAGS) $< $(OBJ) $(LDFLAGS)

# Reciper for file in src
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(SOURCE_DIR)/%.h
	mkdir -p $(dir $@)
	$(CXX) -c -o $@ $(CXXFLAGS) $<

# Recipe for local files
$(OBJ_DIR)/%.o: %.cpp %.h
	mkdir -p $(dir $@)
	$(CXX) -c -o $@ $(CXXFLAGS) $<


.PHONY: clean

clean:
	rm -r $(OBJ) *.exe
