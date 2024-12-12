CXXFLAGS = -O3 --std=c++17 -pedantic\
		 -Wall -Wextra -Werror\
		 -march=native -flto\
		 -fopenmp\
		 -ffast-math -DNDEBUG

AFLAGS = -larmadillo

# Choosing some option for different compilers

ifeq ($(COMPILER), CLANG_LIBCPP)
	CC = clang++-19
	CXXFLAGS+= -stdlib=libc++ -fexperimental-library
else ifeq ($(COMPILER), CLANG_GNUCPP)
	CC = clang++-19
	CXXFLAGS += -stdlib=libstc++
	AFLAGS += -ltbb
else ifeq ($(COMPILER), GNU)
	CC = g++
	AFLAGS += -ltbb
else # Defaults to CLANG_LIBCPP
	CC = clang++-19
	CXXFLAGS+= -stdlib=libc++ -fexperimental-library
endif

ifeq ($(USE_MKL), 1)
	AFLAGS = -DARMA_DONT_USE_WRAPPER  -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
	CXXFLAGS += -I"${MKLROOT}/include"
endif

all: simulator.exe

simulator.exe : system.o NSL_SIMULATOR.o random.o
	$(CC) $(CXXFLAGS) system.o random.o NSL_SIMULATOR.o -o simulator.exe  $(AFLAGS)

NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CXXFLAGS)

system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CXXFLAGS)

# particle.o : particle.cpp particle.h
#	# $(CC) -c particle.cpp -o particle.o $(CXXFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CXXFLAGS)

clean :
	rm *.o simulator.exe

remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
