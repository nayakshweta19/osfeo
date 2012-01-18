# Microsoft Developer Studio Generated NMAKE File, Based on api.dsp
!IF "$(CFG)" == ""
CFG=api - Win32 Debug
!MESSAGE No configuration specified. Defaulting to api - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "api - Win32 Release" && "$(CFG)" != "api - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "api.mak" CFG="api - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "api - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "api - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "api - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\api.lib"


CLEAN :
	-@erase "$(INTDIR)\Channel.obj"
	-@erase "$(INTDIR)\ChannelAddress.obj"
	-@erase "$(INTDIR)\Message.obj"
	-@erase "$(INTDIR)\MovableObject.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\api.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\api.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\api.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\api.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Channel.obj" \
	"$(INTDIR)\ChannelAddress.obj" \
	"$(INTDIR)\Message.obj" \
	"$(INTDIR)\MovableObject.obj"

"$(OUTDIR)\api.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "api - Win32 Debug"

OUTDIR=.\..\..\lib
INTDIR=.\..\..\obj\api
# Begin Custom Macros
OutDir=.\..\..\lib
# End Custom Macros

ALL : "$(OUTDIR)\api.lib"


CLEAN :
	-@erase "$(INTDIR)\Channel.obj"
	-@erase "$(INTDIR)\ChannelAddress.obj"
	-@erase "$(INTDIR)\Message.obj"
	-@erase "$(INTDIR)\MovableObject.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\api.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

"$(INTDIR)" :
    if not exist "$(INTDIR)/$(NULL)" mkdir "$(INTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MLd /W3 /GX /Od /I "..\..\src\api" /I "..\..\src" /I "..\..\src\tagged" /I "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\api.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\api.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"$(OUTDIR)\api.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Channel.obj" \
	"$(INTDIR)\ChannelAddress.obj" \
	"$(INTDIR)\Message.obj" \
	"$(INTDIR)\MovableObject.obj"

"$(OUTDIR)\api.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("api.dep")
!INCLUDE "api.dep"
!ELSE 
!MESSAGE Warning: cannot find "api.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "api - Win32 Release" || "$(CFG)" == "api - Win32 Debug"
SOURCE=..\..\Src\api\elementAPI.cpp

"$(INTDIR)\Channel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Src\api\packages.cpp

"$(INTDIR)\ChannelAddress.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

