# Microsoft Developer Studio Generated NMAKE File, Based on triangle.dsp
!IF "$(CFG)" == ""
CFG=triangle - Win32 Debug
!MESSAGE No configuration specified. Defaulting to triangle - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "triangle - Win32 Release" && "$(CFG)" != "triangle - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "triangle.mak" CFG="triangle - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "triangle - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "triangle - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "triangle - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\triangle.lib"


CLEAN :
	-@erase "$(INTDIR)\Channel.obj"
	-@erase "$(INTDIR)\ChannelAddress.obj"
	-@erase "$(INTDIR)\Message.obj"
	-@erase "$(INTDIR)\MovableObject.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\triangle.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\triangle.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\triangle.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\triangle.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Channel.obj" \
	"$(INTDIR)\ChannelAddress.obj" \
	"$(INTDIR)\Message.obj" \
	"$(INTDIR)\MovableObject.obj"

"$(OUTDIR)\triangle.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "triangle - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\triangle
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\triangle.lib"


CLEAN :
	-@erase "$(INTDIR)\Channel.obj"
	-@erase "$(INTDIR)\ChannelAddress.obj"
	-@erase "$(INTDIR)\Message.obj"
	-@erase "$(INTDIR)\MovableObject.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\triangle.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\actor\objectBroker" /I "..\..\src\actor\address" /I "..\..\src\actor\message" /I "..\..\src\actor\actor" /I "..\..\src\actor\channel" /I "..\..\src\matrix" /I "..\..\src\domain\domain" /I "..\..\src\tagged\storage" /I "..\..\src" /I "..\..\src\tagged" /I "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\actor.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\triangle.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\triangle.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Channel.obj" \
	"$(INTDIR)\ChannelAddress.obj" \
	"$(INTDIR)\Message.obj" \
	"$(INTDIR)\MovableObject.obj"

"$(OUTDIR)\triangle.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("triangle.dep")
!INCLUDE "triangle.dep"
!ELSE 
!MESSAGE Warning: cannot find "triangle.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "triangle - Win32 Release" || "$(CFG)" == "triangle - Win32 Debug"
SOURCE=..\..\Src\actor\channel\Channel.cpp

"$(INTDIR)\Channel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\actor\address\ChannelAddress.cpp

"$(INTDIR)\ChannelAddress.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\actor\message\Message.cpp

"$(INTDIR)\Message.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\actor\actor\MovableObject.cpp

"$(INTDIR)\MovableObject.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

