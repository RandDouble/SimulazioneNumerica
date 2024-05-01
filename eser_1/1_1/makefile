PROGNAME= main_1_1.exe
MAIN= main_1_1.cpp 
SRC_FILES = random.cpp monte_carlo.cpp utilities.cpp

# Compiler and compiler flags
# Default compiler is Clang++-17 with libc++ as standard library implementation
INCFLAGS= -I../../src
LDFLAGS =
CXXFLAGS= $(INCFLAGS)\
		-Wall -Wextra -Werror\
		-std=c++17 -g3\
		-O2
# O2 because it can show more warning due to flow variable analisys, without using too much times
# c++17 because it has some niceties, like in better algorithm and numeric, enabling parallelism

# Choosing some option for different compilers

ifeq ($(COMPILER), CLANG_LIBCPP)
	CXX = clang++-17
	CXXFLAGS += -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CXX = clang++-17
	CXXFLAGS += -stdlib=libstc++
	LDFLAGS += -ltbb
else ifeq ($(COMPILER), GNU)
	CXX = g++
	CXXFLAGS += -Wno-unknown-pragmas
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