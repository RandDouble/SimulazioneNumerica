# @author: Stefano Pilosio
# @scope: Laboratorio di Simulazione Numerica 2023 / 2024

PROGNAME= test_random.exe
MAIN= test_random.cpp
SRC_FILES = old_random.cpp
HEADER_FILES =
LOCAL_SRC_FILE = random_gen_sim.cpp utilities_local.cpp

# Compiler and compiler flags
# Default compiler is Clang++-19 with libc++ as standard library implementation
INCFLAGS= -I../src
SOURCE_DIR= ../src
OBJ_DIR= ./obj

CXXFLAGS= $(INCFLAGS) -Wall -Wextra -O3 -flto -std=c++20 -march=native -DUSE_PARALLEL
# O2 because it can show more warning due to flow variable analisys
# c++19 because it has some niceties, like in better algorithm and numeric, enabling parallelism


ifeq ($(COMPILER), CLANG_LIBCPP) # Only because it support c++20 needed for other projects
	CXX = clang++-19
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CXX = clang++-19
	CXXFLAGS += -stdlib=libstc++
	LDFLAGS += -ltbb
else ifeq ($(COMPILER), GNU)
	CXX = g++
	LDFLAGS += -ltbb
else # Defaults to CLANG_LIBCPP, flags to set library and use parallel algorithms
	CXX = clang++-18
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
endif

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
