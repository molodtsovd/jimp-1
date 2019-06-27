# OMP = -fopenmp
OMP = 
CPPCOMP = g++
CCOMP = gcc
COPT = -c -fPIC -Ofast -march=native -w $(OMP)
#FORTCOMP = gfortran
#FORTOPT = -fopenmp -c -Ofast -march=native -extend-source -ffree-line-length-none

TARGET = tt3d.so
SRCDIR = ./
OBJDIR = ./


.PHONY: all clean

all: $(TARGET)
clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/$(TARGET)


OBJS =  $(OBJDIR)/time_3d.o\
	$(OBJDIR)/time_2d.o\
	$(OBJDIR)/FPT3D.o

## Link objs
$(TARGET): $(OBJS) 
	$(CPPCOMP) $(OMP) -shared -Wl,-soname,$(OBJDIR)/$(TARGET) -o $(OBJDIR)/$(TARGET) $(OBJS)


## Compile

$(OBJDIR)/time_3d.o: $(SRCDIR)/time_3d.c
	$(CCOMP) $(COPT) $(SRCDIR)/time_3d.c -o $(OBJDIR)/time_3d.o

$(OBJDIR)/time_2d.o: $(SRCDIR)/time_2d.c
	$(CCOMP) $(COPT) $(SRCDIR)/time_2d.c -o $(OBJDIR)/time_2d.o

$(OBJDIR)/FPT3D.o: $(SRCDIR)/FPT3D.cpp
	$(CPPCOMP) $(COPT) $(SRCDIR)/FPT3D.cpp -o $(OBJDIR)/FPT3D.o






















