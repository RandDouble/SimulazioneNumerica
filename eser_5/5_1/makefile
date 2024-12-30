# @author: Stefano Pilosio
# @scope: Laboratorio di Simulazione Numerica 2023 / 2024

PROGNAME= main_5_1.exe
MAIN= main_5_1.cpp
SRC_FILES = random.cpp utilities.cpp vector3D.cpp
HEADER_FILES = metropolis.h
LOCAL_SRC_FILE = wave_function.cpp options.cpp
EXTERNAL_HEADER_FILE = json.hpp

# Compiler and compiler flags
# Default compiler is Clang++-19 with libc++ as standard library implementation
INCFLAGS= -I../../src -I../../external_include
CXXFLAGS= $(INCFLAGS) -g -Wall -Wextra -Werror -O3 -std=c++17 -DUSE_PARALLEL # -DNDEBUG
# O2 because it can show more warning due to flow variable analisys
# c++19 because it has some niceties, like in better algorithm and numeric, enabling parallelism

ifeq ($(COMPILER), CLANG_LIBCPP) # Only because it support c++20 needed for other projects
	CXX = clang++-19
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CXX = clang++-19
	CXXFLAGS += -stdlib=libstdc++
	LDFLAGS += -ltbb
else ifeq ($(COMPILER), GNU)
	CXX = g++
	LDFLAGS += -ltbb
else # Defaults to CLANG_LIBCPP, flags to set library and use parallel algorithms
	CXX = clang++-19
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
endif

SOURCE_DIR= ../../src
OBJ_DIR= ./obj

SRC= $(addprefix $(SOURCE_DIR)/, $(SRC_FILES))
HEAD= $(addprefix $(SOURCE_DIR)/, $(HEADER_FILES))

OBJ= $(addprefix $(OBJ_DIR)/,$(addsuffix .o , $(basename $(SRC_FILES))))
SRC += $(LOCAL_SRC_FILE)
OBJ += $(addprefix $(OBJ_DIR)/, $(addsuffix .o , $(basename $(LOCAL_SRC_FILE))))


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


.PHONY: clean

clean:
	rm -r $(OBJ) *.exe
