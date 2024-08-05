#CAPD installation settings
CAPD = ~/CAPD/build/bin
INCLUDE = `$(CAPD)/capd-config --cflags` -I. -I./include/
LIBS =  `$(CAPD)/capd-config --libs` -lpthread
# compiler 
CXX = `$(CAPD)/capd-config --variable=capd_cxx` -O2 $(INCLUDE)

# programs to compile 
CSOURCES = $(wildcard *.cc)
CPPSOURCES = $(wildcard *.cpp)
SOURCES = $(CSOURCES) $(CPPSOURCES)

OBJS = $(patsubst %.cc,obj/%.o,$(CSOURCES)) 
DEPS = $(patsubst %.cpp,dep/%.d,$(CPPSOURCES)) $(patsubst %.cc,dep/%.d,$(CSOURCES))
PROGS = $(patsubst %.cpp,%,$(CPPSOURCES))

all: $(PROGS)

$(PROGS): $(DEPS) $(OBJS) 

%: obj/%.o $(OBJS)
	$(CXX) $< $(OBJS) -o $@ $(LIBS)

include $(DEPS)

dep/%.d: %.c*
	$(CXX) -MM -MT obj/$*.o $< > $@

obj/%.o: %.c* dep/%.d
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o dep/*.d $(PROGS)





