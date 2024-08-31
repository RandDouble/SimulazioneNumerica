# @author: Stefano Pilosio
# @scope: Laboratorio di Simulazione Numerica 2023 / 2024

PROGNAME= main_8_2.exe
MAIN= main_8_2.cpp
SRC_FILES = random.cpp utilities.cpp
HEADER_FILES = metropolis.h
LOCAL_SRC_FILE = ../local_src/wave_function.cpp ../local_src/simulated_annealing.cpp

# Compiler and compiler flags
# Default compiler is Clang++-17 with libc++ as standard library implementation
INCFLAGS= -I../../src -I../local_src
CXXFLAGS= $(INCFLAGS) -std=c++17 -Wall -Wextra -Werror -std=c++17 -O3 -DUSE_PARALLEL -DNDEBUG -ffast-math -march=native
# O2 because it can show more warning due to flow variable analisys
# c++17 because it has some niceties, like in better algorithm and numeric, enabling parallelism

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

SRC= $(addprefix $(SOURCE_DIR)/, $(SRC_FILES))
HEAD= $(addprefix $(SOURCE_DIR)/, $(HEADER_FILES))

OBJ= $(addprefix $(OBJ_DIR)/,$(SRC_FILES:.cpp=.o))
SRC += $(LOCAL_SRC_FILE)
OBJ += $(addprefix $(OBJ_DIR)/, $(LOCAL_SRC_FILE:.cpp=.o))


.PHONY: all
all : $(PROGNAME)

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

$(OBJ_DIR)/%.o: ../local_src/%.cpp ../local_src/%.h
	mkdir -p $(dir $@)
	$(CXX) -c -o $@ $(CXXFLAGS) $<

test.exe: test_rng.cpp $(OBJ) $(HEAD)
	$(CXX) -o test.exe $(CXXFLAGS) test_rng.cpp $(OBJ) $(LDFLAGS)


.PHONY: clean

clean:
	rm -r $(OBJ) *.exe
