include ../../Makefile.def

OBJS  = ExperimentalElement.o \
        EEBeamColumn2d.o \
        EEBeamColumn3d.o \
        EEBearing2d.o \
        EEBearing3d.o \
        EEGeneric.o \
        EEInvertedVBrace2d.o \
        EETruss.o \
        EETrussCorot.o \
        EETwoNodeLink.o

all:         $(OBJS)

tcl: 	TclEEBeamColumnCommand.o \
    TclEEBearingCommand.o \
	TclEEGenericCommand.o \
	TclEEInvertedVBraceCommand.o \
	TclEETrussCommand.o \
	TclEETwoNodeLinkCommand.o \
	TclExpElementCommands.o 



# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core example main.o

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
