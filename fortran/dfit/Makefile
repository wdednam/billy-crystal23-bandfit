#------------------------------------------------------------------#
# Makefile for dfit program.                                       #
#                                                                  #
# MDT 1.2005                                                       #
#------------------------------------------------------------------#

# Define fortran compiler here
F90 = ifort -O3 -xCORE-AVX2
#F90 = gfortran -O3
LDF90        = $(F90)

# Basic stuff
UTILEXE = dfit
SRCDIR  = $(PWD)
BINDIR  = $(PWD)
OBJDIR  = $(PWD)
.SILENT:

start: default
.PHONY: start

# Main target
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# Object-file list
OBJ = \
$(OBJDIR)/dfit.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/errvrs.o \
$(OBJDIR)/lesolv.o \
$(OBJDIR)/mxmb.o
# End object-file list

# Module dependencies
DFIT_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/dfit.o : $(SRCDIR)/dfit.f90 $(DFIT_DEP)
dfit.o_name = dfit

$(OBJDIR)/dsp.o : $(SRCDIR)/dsp.f90
dsp.o_name = dsp

$(OBJDIR)/errvrs.o : $(SRCDIR)/errvrs.f90
errvrs.o_name = errvrs

LESOLV_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/lesolv.o : $(SRCDIR)/lesolv.f90 $(LESOLV_DEP)
lesolv.o_name = lesolv

MXMB_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/mxmb.o : $(SRCDIR)/mxmb.f90 $(MXMB_DEP)
mxmb.o_name = mxmb
# End module dependencies

# How to make the objects
$(OBJ):
	$(F90) -c $(SRCDIR)/$($(@F)_name).f90
$(BINDIR)/$(UTILEXE) : $(OBJ)
	$(LDF90) -o $(@) $(OBJ) ; rm -f *.o *.mod >& /dev/null
