# OMP = -fopenmp
OMP = 
CCOMP = gcc
COPT = -c -fPIC -Ofast -w $(OMP) -march=native
FORTCOMP = gfortran
FORTOPT = $(OMP) -c -fPIC -Ofast -ffree-line-length-none #-march=native

TARGET = compro.so
SRCDIR = ./
OBJDIR = ./


.PHONY: all clean

all: $(TARGET)
clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/$(TARGET)


OBJS =  $(OBJDIR)/compro.o\
	$(OBJDIR)/sub_compro.o

## Link objs
$(TARGET): $(OBJS) 
	$(FORTCOMP) $(OMP) -shared -Wl,-soname,$(OBJDIR)/$(TARGET) -o $(OBJDIR)/$(TARGET) $(OBJS)


## Compile

$(OBJDIR)/sub_compro.o: $(SRCDIR)/sub_compro.f90
	$(FORTCOMP) $(FORTOPT) $(SRCDIR)/sub_compro.f90 -o $(OBJDIR)/sub_compro.o

$(OBJDIR)/compro.o: $(SRCDIR)/compro.c
	$(CCOMP) $(COPT) $(SRCDIR)/compro.c -o $(OBJDIR)/compro.o




# 
# CCOMP = gcc
# COPT = -c -Ofast -march=native -w -fopenmp
# FORTCOMP = gfortran
# FORTOPT = -fopenmp -c -Ofast -march=native -extend-source -ffree-line-length-none
# 
# TARGET = compro
# SRCDIR = ./
# OBJDIR = ./
# 
# 
# .PHONY: all clean
# 
# all: $(TARGET)
# clean:
# 	rm $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/$(TARGET)
# 
# 
# OBJS =  $(OBJDIR)/compro.o\
# 	$(OBJDIR)/sub_compro.o
# 
# ## Link objs
# $(TARGET): $(OBJS) 
# 	$(FORTCOMP) -fopenmp -o $(OBJDIR)/$(TARGET) $(OBJS)
# 
# 
# ## Compile
# 
# $(OBJDIR)/sub_compro.o: $(SRCDIR)/sub_compro.f90
# 	$(FORTCOMP) $(FORTOPT) $(SRCDIR)/sub_compro.f90 -o $(OBJDIR)/sub_compro.o
# 
# $(OBJDIR)/compro.o: $(SRCDIR)/compro.c
# 	$(CCOMP) $(COPT) $(SRCDIR)/compro.c -o $(OBJDIR)/compro.o
# 
# 
# 















