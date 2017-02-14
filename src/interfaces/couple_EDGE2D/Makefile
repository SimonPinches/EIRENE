################################
# Detect platform if undefined #
################################
ifeq ($(flavor PLATFORM),undefined)
   PLATFORM=$(shell uname -i)
   ifeq ("$(PLATFORM)","i386")
      PLATFORM=i686
   endif
   $(info PLATFORM detected as $(PLATFORM))
endif
export PLATFORM

###########################################
# Set up compilation specific environment #
###########################################

# set SHELL to bash
# -----------------
SHELL:=/bin/bash

# Set Fortran and C compiler
# --------------------------
F77=pgf90
CC=/usr/bin/gcc

# Set compilation flags
# ---------------------
ifeq ("$(DBG)","")
#+PJK
   ifeq ($(PLATFORM),x86_64)
      FFLAGS= -c -O2 -m64 -mcmodel=medium -tp=sandybridge
   else
      FFLAGS= -c -O2
   endif
else
   ifeq ($(PLATFORM),x86_64)
      FFLAGS= -c -g -m64 -mcmodel=medium -tp=sandybridge
   else
      FFLAGS= -c -g
   endif
#-PJK
endif
CFLAGS= -c -O2

# Set compilation options
# -----------------------
OPTS= -Mextend
#OPTS= -Mextend -Mnomain
#OPTS= -Mnolist -Mdepchk -Minform=inform -Minfo -Mchkptr -Mchkfpstk -Mfptrap -Mbounds
#OPTS= -Mnolist -Mdepchk -Minform=inform -Minfo -Mchkfpstk -Mfptrap -Mbounds
COPTS=

# Set Libraries to link with
# --------------------------
#LIBS= -L/usr/pgi/linux86/lib-glibc-211 -L/usr/local/lib -l gr -l gks  -L/usr/X11R6/lib/ -lX11 -lXt
#LIBS= -L/usr/pgi/linux86/lib-glibc-211 -L/usr/local/lib  -L/usr/X11R6/lib/  -lX11 -lXt
LIBS =

# Set include directory
# ---------------------
IDIR=

# set path to search for Prerequisites
# ------------------------------------
VPATH=./src:./modules:./

############################################
# Set up PREPROCESSOR specific environment #
############################################

# get GIT release tag
# -------------------
EIR_PREPROC_GITTAG:=$(shell git describe --tags --abbrev=0 --match Release-v[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]])

# set Release Date string (eg. 01-JAN-01)
# ---------------------------------------
EIR_PREPROC_RELEASEDATE:="'$(shell a=$(EIR_PREPROC_GITTAG); date --date=$${a:13:2}$${a:11:2}$${a:9:2} +"%d-%b-%y" | tr "[:lower:]" "[:upper:]" )'"

# set wanted JETTO version
# ------------------------
EIR_PREPROC_WANTEDJETTOVERSION:="'$(shell a=$(EIR_PREPROC_GITTAG); echo j$${a:9:6})'"

# set actual GIT branch name
# --------------------------
dummy:=$(shell git branch --no-color | grep [*])
EIR_PREPROC_GITBRANCHNAME:="'$(shell a="$(dummy)"; echo $${a:2})'"

# set actual GIT status string
# ----------------------------
ifeq ($(shell git status --porcelain),)
EIR_PREPROC_GITSTATUS:="'clean'"
else
EIR_PREPROC_GITSTATUS:="'modyf'"
endif

# set actual GIT commit info
# --------------------------
EIR_PREPROC_GITSHAKEY:="'$(shell git --git-dir=./.git rev-parse --verify HEAD)'"

# set actual GIT repository
# -------------------------
EIR_PREPROC_GITREPOSITORY:="'$(shell pwd)'"

# set hostname used for compilation
# ---------------------------------
EIR_PREPROC_HOSTNAME:="'$(shell hostname)'"

# format actual GIT tag string
# ----------------------------
EIR_PREPROC_GITTAG:="'$(EIR_PREPROC_GITTAG)'"

######################################
# Set up EIRENE specific environment #
######################################
  
# Get build name
# --------------
ESRCPATH:=$(shell pwd)
EBUILDNAME:=$(notdir $(ESRCPATH))
$(info )
$(info Running Makefile for EIRENE build: $(EBUILDNAME))
$(info )

# Get last EIRENE git release tag of format Release-vDDMMYY
# ---------------------------------------------------------
EGITTAG:=$(shell git describe --tags --abbrev=0 --match Release-v[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]] )
$(info Last EIRENE git release TAG: $(EGITTAG))

# Get EIRENE version identifier from last git release tag
# -------------------------------------------------------
EGITTAGVER:=$(strip $(subst Release-, ,${EGITTAG}))
$(info EIRENE source is based on release version: ${EGITTAGVER})
$(info Building EIRENE library in: ${ESRCPATH})
$(info )

# Check if this is a release build under $GUSR
# --------------------------------------------
WHOAMI=$(shell whoami)
ifeq ($(GUSR), $(WHOAMI))
  ifeq ($(EGITTAGVER), $(EBUILDNAME))
    ifneq ($(JINTRACRELEASE), $(EGITTAG))
      $(error Environment varibale JINTRACRELEASE is not set to $(EGITTAG)!)
    endif
  endif
endif

# Check if this is a bugfix build under $GUSR and adjust EIR_PREPROC_GITREPOSITORY
# --------------------------------------------------------------------------------
ifeq ($(GUSR), $(WHOAMI))
  RELEASEVERSION:=$(shell b=$(EBUILDNAME); a=$${b:0:8}; if [[ $$a =~ v[[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]][[:digit:]]_ ]]; then echo "true"; else echo "false"; fi)
  ifeq ($(RELEASEVERSION),true)
    EBUGFIX:=$(shell a=$(EBUILDNAME); echo $${a:8:6})
    ifeq ($(EBUGFIX),bugfix)
      EIR_PREPROC_GITREPOSITORY:="'$(dir $(ESRCPATH))$(shell a=$(EBUILDNAME); echo $${a:0:7})'"
      $(info Compiling bugfix version under user GUSR)
      $(info Adjusted EIR_PREPROC_GITREPOSITORY to $(EIR_PREPROC_GITREPOSITORY))
      $(info )
    endif
  endif
endif

# program name
# ------------
#+PJK
ifeq ($(PLATFORM),x86_64)
   PROG=libeirene_64.a
else
   PROG=libeirene.a
endif
#-PJK
# Executable target name
# ----------------------
TARG=eirene

# Delete default suffixes
# -----------------------
.SUFFIXES:

# Set suffixes
# ------------
.SUFFIXES: .o .f .f90

# Set PHONY targets
# -----------------
.PHONY: all clean realclean force_git_status

# Set module list
# ---------------
MODULES = precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o mod_eckstein.o

# Set object list
# ---------------
OBJECTS= photons_dummy.o main-routines.o startup-routines.o particle-tracing.o volume-processes.o surface-processes.o file-handling.o sampling.o scoring.o output.o geometry.o time-routines.o plotting.o iterate.o tetrahedrons.o assistant.o mathematics.o diagno.o broadcast_dummy.o couple_Dummy.o nag-dummy.o grsoft-dummy.o xdr-dummy.o user_Edge2d.o git_Edge2d_info.o


# make 'all' option
# -----------------
all: $(PROG)

# Rule to build eirene library
# ----------------------------
$(PROG): $(MODULES) $(OBJECTS)
	ar vr $(PROG) $?
#	ar vr $(PROG) $(MODULES) $(OBJECTS)

# Rule to build executable
# ------------------------
#$(PROG): $(MODULES) $(OBJECTS)
#	$(F77) -o $(TARG) $(MODULES) $(OBJECTS) $(LIBS)

# Standard rules for source files
# -------------------------------
.c.o:
	$(CC) $(CFLAGS) $(COPTS) $?
.f.o:
	$(F77) $(FFLAGS) $(OPTS) $?
.f90.o:
	$(F77) $(FFLAGS) $(OPTS) $?

# set rule for git_Edge2d_info.f (include c-preprocessor to pass GIT information to executable)
# ---------------------------------------------------------------------------------------------
git_Edge2d_info.o: git_Edge2d_info.f git_status_old
	$(F77) $(FFLAGS) $(OPTS) -Mextend \
	-D__WANTEDJETTOVERSION=$(EIR_PREPROC_WANTEDJETTOVERSION) \
	-D__GITRELEASEDATE=$(EIR_PREPROC_RELEASEDATE) \
	-D__GITRELEASETAG=$(EIR_PREPROC_GITTAG) \
	-D__GITBRANCHNAME=$(EIR_PREPROC_GITBRANCHNAME) \
	-D__GITSTATUSSTRING=$(EIR_PREPROC_GITSTATUS) \
	-D__GITSHAKEY=$(EIR_PREPROC_GITSHAKEY) \
	-D__GITREPOSITORY=$(EIR_PREPROC_GITREPOSITORY) \
	-D__HOSTNAME=$(EIR_PREPROC_HOSTNAME) \
	-Mpreprocess $<

# rule to check GIT status
# ------------------------
git_status_old: force_git_status
	@if [ ! -f git_status_old ] ; then \
	    touch git_status_old; \
	fi
	@echo $(EIR_PREPROC_WANTEDJETTOVERSION) $(EIR_PREPROC_RELEASEDATE) \
	      $(EIR_PREPROC_GITTAG) $(EIR_PREPROC_GITBRANCHNAME) $(EIR_PREPROC_GITSTATUS) \
	      $(EIR_PREPROC_GITSHAKEY) $(EIR_PREPROC_GITREPOSITORY) $(EIR_PREPROC_HOSTNAME) > git_status_new
	@a=$$(diff git_status_new git_status_old); \
	if [ $$? != 0 ]; then mv git_status_new git_status_old; echo "Repository status has changed!"; fi

# make 'realclean' option (remove objects, modules and libraries)
# ---------------------------------------------------------------
realclean:
	rm -f *.o *.mod *.a

# make 'clean' option (just remove object files)
# Currently remove also *.mod files for safe dual build of 32 and 64 bit
# ----------------------------------------------------------------------
clean:
	rm -f *.o *.mod

# Dependencies
# ------------
precision.o:
parmmod.o: precision.o
photons_dummy.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
main-routines.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
startup-routines.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
particle-tracing.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
volume-processes.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
surface-processes.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o  mod_eckstein.o
file-handling.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
sampling.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
scoring.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
output.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
geometry.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
time-routines.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
plotting.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
iterate.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
tetrahedrons.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
assistant.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
mathematics.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
diagno.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
user_C-Mod.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
user.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
user_Edge2d.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
broadcast_dummy.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
couple_Dummy.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
nag-dummy.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
grsoft-dummy.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
xdr-dummy.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
outusr.o: precision.o parmmod.o braeir.o braspoi.o cadgeo.o cai.o ccona.o ccoupl.o cestim.o cfplk.o cgeom.o cgrid.o cgrptl.o cinit.o clast.o clgin.o clmsur.o clogau.o comnnl.o comprt.o comsig.o comsou.o comspl.o comusr.o comxs.o coutau.o cpes.o cpl3d.o cplmsk.o cplot.o cpolyg.o crand.o crech.o cref.o crefmod.o csdvi_bgk.o csdvi_cop.o csdvi.o cspei.o cspez.o cstep.o ctetra.o ctext.o ctrcei.o ctrig.o ctsurf.o cupd.o czt1.o eirbra.o module_avltree.o
eirene_xs_util.o: mod_eirene_xs.o
