#Boris Jeremic 06Mar2009

DEBUG_MODE = DEBUG
#DEBUG_MODE = NO_DEBUG

#Compiler = GNU
Compiler = PGI

BASE  = /usr
#  PUT YOUR HOME DIRECTOREY HERE

PWD := $(shell pwd)

  PATH_TO_NewTemplate3Dep       = $(PWD)

include ../../Makefile.Compilers



MACHINE_INCLUDES        = -I/usr/include 
# this file contains all the OpenSees1/ModifiedOpenSeesServices includes


FE_INCLUDES  = -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/matrix \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/tagged \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/material \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/material/nD \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/recorder \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/handler \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/actor \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/analysis \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/ \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/actor/channel \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/actor/actor \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/element \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/domain \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/domain/domain \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/domain/partitioner \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/domain/component \
               -I$(PATH_TO_NewTemplate3Dep)/../../ModifiedOpenSeesServices/actor/objectBroker \
               -I$(PATH_TO_NewTemplate3Dep)/../../nDarray \
               -I$(PATH_TO_NewTemplate3Dep)/. \



INCLUDES =  $(MACHINE_INCLUDES)  $(FE_INCLUDES)


# COMPILATION

.SUFFIXES:
.SUFFIXES:	.C .c .f .f90 .cpp .o .cpp

#
# %------------------%
# | Default command. |
# %------------------%
#
.DEFAULT:
	@$(ECHO) "Unknown target $@, try:  make help"
#
# %-----------------------------------------------%
# |  Command to build .o files from source files. |
# %-----------------------------------------------%
#

.cpp.o: 
#	@$(ECHO) Making $@ from $< $@ with $(CC++) $(C++FLAGS) $(INCLUDES) -c $<
	$(CC++) $(C++FLAGS) $(INCLUDES) -c $<


#include ../../Makefile.def


OBJS =  MaterialParameter.o \
        ElasticState.o \
        Isotropic_Elastic.o \
        elnp_Elastic.o \
        PressureDependent_Elastic.o \
        DM04_Elastic.o \
        YieldFunction.o \
        VM_YF.o \
        RMC_YF.o \
        DP_YF.o \
        CC_YF.o \
        DM04_YF.o \
        PlasticFlow.o \
        VM_PF.o \
        RMC_PF.o \
        DP_PF.o \
        CC_PF.o \
        DM04_PF.o \
        ScalarEvolution.o \
        Linear_Eeq.o \
        CC_Ev.o \
        TensorEvolution.o \
        Linear_Eij.o \
        AF_Eij.o \
        DM04_alpha_Eij.o \
        DM04_z_Eij.o \
        NewTemplate3Dep.o \
        SANISAND_YF.o \
        SANISAND_Elastic.o \
        SANISAND_alpha_Eij.o \
        SANISAND_PF.o \
        SANISAND_p0_bar.o \
        SANISAND_z_Eij.o \


#        TclNewTemplate3Dep.o \

LIBRARY_NewTemplate_3D_EP  = ../../lib/libnewtemplate3dep.a


################### TARGETS ########################
all: $(OBJS) 
	$(AR) $(ARFLAGS) $(LIBRARY_NewTemplate_3D_EP) $(OBJS)

#libmatmod.a: $(OBJS)
#	$(AR) $(ARFLAGS) $(Template_3D_EP) $(OBJS)




#3D_El_Pl        = 3D_EL_PL
#
#Constitutive_Testing:  Constitutive_Testing.o 
#	$(LINKER) $(LINKFLAGS) Constitutive_Testing.o \
#	$(INTERPRETER_LIBS_TCL) \
#	$(FE_LIBRARY) \
#	$(MACHINE_LINKLIBS) \
#	$(MACHINE_NUMERICAL_LIBS) \
#	$(MACHINE_SPECIFIC_LIBS) \
#	-o Constitutive_Testing
#
#3D_El_Pl:  Constitutive_Testing.o 
#	$(LINKER) $(LINKFLAGS) Constitutive_Testing.o \
#	$(INTERPRETER_LIBS_TCL) \
#	$(FE_LIBRARY) \
#	$(MACHINE_LINKLIBS) \
#	$(MACHINE_NUMERICAL_LIBS) \
#	$(MACHINE_SPECIFIC_LIBS) \
#	-o $(3D_El_Pl)
#
#Bardet:  BardetConstitutiveTesting.o 
#	$(LINKER) $(LINKFLAGS) BardetConstitutiveTesting.o ../BardetConstraint.o \
#	$(FE_LIBRARY) \
#	$(MACHINE_LINKLIBS) \
#	$(MACHINE_NUMERICAL_LIBS) \
#	$(MACHINE_SPECIFIC_LIBS) \
#	-o BardetConstitutiveTesting
#
##	$(NDARRAY_LIBRARY) \
#
#
#
#FT3Dep_test: $(OBJS) FT3Dep_test.o
#	$(LINKER) $(LINKFLAGS) FT3Dep_test.o $(OBJS) \
#	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY)\
#	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
#	-o FT3Dep_test


clean:
	$(RM) $(RMFLAGS)  *.o *.bak

cleanall:
	$(RM) $(RMFLAGS)  *.o *.bak
	$(RM) $(RMFLAGS) $(LIBRARY_NewTemplate_3D_EP)



# DO NOT DELETE THIS LINE -- make depend depends on it.




