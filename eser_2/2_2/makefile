PROGNAME= main_2_2.exe
MAIN= main_2_2.cpp 
SRC_FILES = random.cpp utilities.cpp random_walk.cpp vector3D.cpp
# LOCAL_SRC_FILE = buffon.cpp
# Compiler and compiler flags
CXX=clang++-17 # Only because it support c++20 needed for other projects
INCFLAGS= -I../../src
CXXFLAGS= $(INCFLAGS) -Wall -Wextra -Werror -O2 -std=c++17
# O2 because it can show more warning due to flow variable analisys
# c++17 because it has some niceties, like in better algorithm and numeric, enabling parallelism

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
	$(CXX) -o $@ $(CXXFLAGS) $^ 

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