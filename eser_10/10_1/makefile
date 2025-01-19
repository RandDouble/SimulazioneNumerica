# @author: Stefano Pilosio
# @scope: Laboratorio di Simulazione Numerica 2023 / 2024

PROGNAME= main_10_1.exe
MAIN= main_10_1.cpp
SRC_FILES = random.cpp utilities.cpp
HEADER_FILES = initializer.h mutation.h population.h swap_functions.h
LOCAL_SRC_FILE =
LOCAL_HEADER_FILE =
EXTERNAL_HEADER_FILE = json.hpp

# Compiler and compiler flags
CXX = mpicxx # Change this line to change compiler

CPPFLAGS += -I../../src -I../../external_include
CXXFLAGS += -Wall -Wextra -g -O3 -std=c++17 -pedantic -march=native -DNDEBUG
LDFLAGS += -larmadillo

# O2 because it can show more warning due to flow variable analisys
# c++19 because it has some niceties, like in better algorithm and numeric, enabling parallelism

# Disabled control structure, in this case I want to use only MPICH
# Changed to OpenMPI due to compiler bug in Ubuntu 24.04 LTS WSL2

# ifeq ($(COMPILER), CLANG_LIBCPP) # Only because it support c++20 needed for other projects
# 	CXX = clang++-19
# 	CXXFLAGS += -stdlib=libc++ -fexperimental-library
# else ifeq ($(COMPILER), CLANG_GNUCPP)
# 	CXX = clang++-19
# 	LDFLAGS += -ltbb
# else ifeq ($(COMPILER), GNU)
# 	CXX = g++
# 	LDFLAGS += -ltbb
# else ifeq ($(COMPILER), MPICH)
# 	CXX = mpicxx
# else # Defaults to CLANG_LIBCPP, flags to set library and use parallel algorithms
# 	CXX = clang++-19
# 	CXXFLAGS += -stdlib=libc++ -fexperimental-library
# endif



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

$(PROGNAME): $(MAIN) $(OBJ) $(HEAD)
	$(CXX) $(CPPFLAGS) -o $@ $(CXXFLAGS) $< $(OBJ) $(LDFLAGS)

# Reciper for file in src
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(SOURCE_DIR)/%.h
	mkdir -p $(dir $@)
	$(CXX) -c $(CPPFLAGS) -o $@ $(CXXFLAGS) $<

# Recipe for local files
$(OBJ_DIR)/%.o: %.cpp %.h
	mkdir -p $(dir $@)
	$(CXX) -c $(CPPFLAGS) -o $@ $(CXXFLAGS) $<


.PHONY: clean
clean:
	rm -r $(OBJ) *.exe

.PHONY: exec

exec: all
	mpiexec -n 6 ./$(PROGNAME)
	mkdir -p ../results/multi_thread/
	mv ../output/* ../results/multi_thread/
# mv ../results/best_route.csv ../results/multi_thread/

.PHONY: exec_single
exec_single: all
	mpiexec -n 1 ./$(PROGNAME)
	mkdir -p ../results/single_thread
	mv ../output/* ../results/single_thread/

# mv ../results/best_route.csv ../results/single_thread/best_route_single_thread.csv
