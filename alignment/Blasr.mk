PBINCROOT ?= /srv/depot/jdrake/p4/bioinformatics/lib/cpp
HDFINC ?= /srv/depot/jdrake/p4/assembly/seymour/dist/common/include

INCDIRS = -I$(PBINCROOT)/alignment \
		  -I$(PBINCROOT)/pbdata \
		  -I$(PBINCROOT)/hdf \
		  -I$(HDFINC) \
	      -I../common 

CXXFLAGS := -std=c++0x -Wall -Wuninitialized -Wno-div-by-zero \
			-pedantic -c -fmessage-length=0 -MMD -MP

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)
LIBS := -lblasr -lpbdata -lpbihdf -lpthread 

all : OPTIMIZE = -O3

debug : OPTIMIZE = -g

profile : OPTIMIZE = -Os -pg

all debug profile: blasr

blasr: Blasr.o
	$(CXX) $(LIBDIRS) -static -o $@ $^ $(LIBS)

Blasr.o: Blasr.cpp
	$(CXX) $(CXXFLAGS) $(OPTIMIZE) $(INCDIRS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -o $@ $<

.INTERMEDIATE: $(OBJS)

clean: 
	rm -f blasr
	rm -f $(OBJS) $(DEPS)

-include $(DEPS)
