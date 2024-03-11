PROGNAME= main_1_1.exe
MAIN= main_1_1.cpp 
SRC_FILES = random.cpp monte_carlo.cpp utilities.cpp

# Compiler and compiler flags
CXX=clang++-17 # Only because it support c++20 needed for other projects
INCFLAGS= -I../../src
CXXFLAGS= $(INCFLAGS)\
 		-Wall\
  		-Wextra\
   		-Werror\
	 	-std=c++17 -g3\
    	-O2\
# O2 because it can show more warning due to flow variable analisys
# c++17 because it has some niceties, like in better algorithm and numeric, enabling parallelism

SOURCE_DIR= ../../src
OBJ_DIR= ./obj
SRC= $(addprefix $(SOURCE_DIR), $(SRC_FILES))
OBJ= $(addprefix $(OBJ_DIR)/,$(addsuffix .o , $(basename $(SRC_FILES))))

$(PROGNAME): $(MAIN) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $^ 

$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) -c -o $@ $(CXXFLAGS) $<

.PHONY: clean

clean:
	rm -r $(OBJ) *.exe