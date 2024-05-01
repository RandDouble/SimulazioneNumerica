DEPS:= eser_1/1_1/main_1_1.exe\
	eser_1/1_2/main_1_2.exe\
	eser_1/1_3/main_1_3.exe\
	eser_2/2_1/main_2_1.exe\
	eser_2/2_2/main_2_2.exe\
	eser_3/3_1/main_3_1.exe\
	NSL_SIMULATOR/SOURCE/simulator.exe\
	eser_5/5_1/main_5_1.exe

# Set COMPILER to choose a compiler platform,
# Possible values are:
# - CLANG_LIBCPP to compiler with clang-17 and libc++
# - CLANG_GNUCPP to compiler with clang-17 and listdbc++
# - GNU to compiler with g++

COMPILER=GNU
export COMPILER

.PHONY: all

all: $(DEPS)

%.exe:
	$(MAKE) -C $(dir $@) -j

.PHONY: clean

clean:
	$(foreach var, $(DEPS), $(MAKE) -C $(dir $(var)) clean;)
