PROGNAME= main_2_2.exe
MAIN= main_2_2.cpp 
SRC_FILES = random.cpp utilities.cpp random_walk.cpp vector3D.cpp
# LOCAL_SRC_FILE = buffon.cpp
# Compiler and compiler flags
INCFLAGS= -I../../src
LDFLAGS=
CXXFLAGS= $(INCFLAGS) -Wall -Wextra -Werror -O3 -std=c++17 -ffast-math -march=native -fopenmp -DUSE_PARALLEL

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
HEAD= $(addprefix $(SOURCE_DIR)/, $(HEADER_FILES))


OBJ= $(addprefix $(OBJ_DIR)/,$(addsuffix .o , $(basename $(SRC_FILES))))
SRC += $(LOCAL_SRC_FILE)
OBJ += $(addprefix $(OBJ_DIR)/, $(addsuffix .o , $(basename $(LOCAL_SRC_FILE))))

.PHONY: all
all : $(PROGNAME)

$(PROGNAME): $(MAIN) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

# Recipe for files in src
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
