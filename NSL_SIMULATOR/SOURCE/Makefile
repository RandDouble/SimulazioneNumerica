CFLAGS = -O3 --std=c++17\
		 -Wall -Wextra -Werror\
		 -march=native\
		 -fopenmp\
		 -ffast-math


AFLAGS = -larmadillo

# Choosing some option for different compilers

ifeq ($(COMPILER), CLANG_LIBCPP)
	CC = clang++-17
	CFLAGS+= -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CC = clang++-17
	CFLAGS += -stdlib=libstc++
	AFLAGS += -ltbb
else ifeq ($(COMPILER), GNU)
	CC = g++
	AFLAGS += -ltbb
else # Defaults to CLANG_LIBCPP
	CC = clang++-17
	CFLAGS+= -stdlib=libc++ -fexperimental-library
endif

all: simulator.exe

simulator.exe : system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) $(CFLAGS) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe  $(AFLAGS)

NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)

system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)

particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o simulator.exe

remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
