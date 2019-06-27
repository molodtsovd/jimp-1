# OMP = -fopenmp
OMP = 
CPPCOMP = g++
COPT = -c -fPIC -Ofast -march=native -w $(OMP)

TARGET = grav.so
SRCDIR = ./
OBJDIR = ./


.PHONY: all clean

all: $(TARGET)
clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/$(TARGET)


OBJS =  $(OBJDIR)/grav.o

## Link objs
$(TARGET): $(OBJS) 
	$(CPPCOMP) $(OMP) -shared -Wl,-soname,$(OBJDIR)/$(TARGET) -o $(OBJDIR)/$(TARGET) $(OBJS)


## Compile

$(OBJDIR)/grav.o: $(SRCDIR)/grav.cpp
	$(CPPCOMP) $(COPT) $(SRCDIR)/grav.cpp -o $(OBJDIR)/grav.o






















