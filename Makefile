#
#
## --------- BEGIN EDIT  --------- ##

## compilers:
CC = g++
MPICXX = mpicxx

VERSION=`cat VERSION`

## plateform options:
## define (uncomment) what suites your needs (nothing needed for linux)
#MAC = 1
#WIN = 1
#MPI = 1

## GSL is always use, check the path below
GSL=1

## uncomment the following to compile the debugging version:
#DEBUG = 1

##install paths:
BIN_INSTALL = /bin/

##library paths:
GSL_PATH = /usr/local/

# you'll have to customize the path to the SPRNG library if you use MPI
#SPRNG_PATH = ../../sprng2.0/

## extra preprocessor macro to remove the runtime generation counter:
#C_OPTS = -DLOW_VERBOSE

## ---------- END EDIT ---------- ##

ifdef DEBUG
BIN_NAME = nemo$(VERSION)D
else
BIN_NAME = nemo$(VERSION)
BIN_NAME = nemo$(VERSION)-locked
endif

#warnings:
W_OPTS = -w

#optimizations:
O_OPTS = -O3 -funroll-loops

#libraries and customizations:

ifdef MPI
    SPRNG = 1
    C_OPTS += -DUSE_MPI
    LD_OPTS = -lmpich
    CCINCL = -I/path/to/mpich/include
endif

ifdef DEBUG
    O_OPTS = -O0
    C_OPTS += -D_DEBUG_ -g #-pg
    W_OPTS = -Wall
endif

ifdef GSL
    C_OPTS += -DHAS_GSL
    CCINCL = -I$(GSL_PATH)include
ifdef WIN
    LD_OPTS = -lstdc++ -static-libgcc -static -lgsl -lgslcblas  #static linking doesn't seem to work anymore
else
ifdef STATIC
    LD_OPTS +=$(GSL_PATH)lib/libgsl.a $(GSL_PATH)lib/libgslcblas.a
else
	LD_OPTS +=-L$(GSL_PATH)lib -lgsl -lgslcblas
endif
endif
endif

ifdef SPRNG
    LD_OPTS += -L$(SPRNG_PATH)lib -lsprng
    CCINCL += -I$(SPRNG_PATH)include
endif

#Platform specific stuff
ifdef WIN
    C_OPTS += -D_WINDOWS_
endif

ifdef MAC
#    O_OPTS += -mtune=G5
#    C_OPTS += -arch i386
endif

ifdef MPI
    BIN_NAME = nemo$(VERSION)_mpi
    LIB_NAME = libnemo$(VERSION)_mpi
    W_OPTS = -Wbrief
    LIB_EXT = .so
    C_OPTS += -DUSE_MPI -DLOW_VERBOSE
    CC = $(MPICXX)
endif

BIN = bin/$(BIN_NAME)
LIB = lib/$(LIB_NAME)$(LIB_EXT)
CCFLAGS = $(C_OPTS) $(W_OPTS) $(O_OPTS)

SRC_PATH = src/
SOURCES = $(shell ls src/*.cc)

OBJ_PATH = src/
OBJECTS = $(shell for file in $(SOURCES);\
                do echo $$file | sed -e "s/\(.*\)\.cc/\1\.o/";\
                done)

all : bin

depend:
	/usr/X11R6/bin/makedepend -Y $(SOURCES)

install: clean all
	cp $(BIN) $(BIN_INSTALL)$(BIN_NAME)
	ln -s $(BIN_INSTALL)$(BIN_NAME) $(BIN_INSTALL)nemo

bin : $(OBJECTS)
	$(CC) -o $(BIN) $(OBJECTS) $(LD_OPTS)

clean :
	rm -f $(OBJECTS); rm -f $(OBJ_PATH)main.o

cleaner: clean
	rm -f $(BIN) $(LIB) *.bak *~ .DS*

$(OBJ_PATH)%.o: $(SRC_PATH)%.cc
	$(CC) $(CCFLAGS) -c $(SRC_PATH)$*.cc -o $(OBJ_PATH)$*.o $(CCINCL)

# DO NOT DELETE THIS LINE -- make depend depends on it.

src/LCEbreed.o: src/output.h src/LCEbreed.h src/lifecycleevent.h src/param.h
src/LCEbreed.o: src/handler.h src/tmatrix.h src/types.h src/simcomponent.h
src/LCEbreed.o: src/fileservices.h src/service.h src/statservices.h
src/LCEbreed.o: src/statrecorder.h src/updaterservices.h
src/LCEbreed.o: src/binarystoragebuffer.h src/metapop.h src/indfactory.h
src/LCEbreed.o: src/individual.h src/ttrait.h src/ttrait_with_map.h
src/LCEbreed.o: src/MPStatHandler.h src/stathandler.h src/binarydataloader.h
src/LCEbreed.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/LCEbreed.o: src/Uniform.h src/simenv.h src/simulation.h
src/LCEbreed.o: src/basicsimulation.h
src/LCEcomposite.o: src/LCEcomposite.h src/LCEbreed.h src/lifecycleevent.h
src/LCEcomposite.o: src/param.h src/handler.h src/tmatrix.h src/output.h
src/LCEcomposite.o: src/types.h src/simcomponent.h src/fileservices.h
src/LCEcomposite.o: src/service.h src/statservices.h src/statrecorder.h
src/LCEcomposite.o: src/updaterservices.h src/binarystoragebuffer.h
src/LCEcomposite.o: src/metapop.h src/indfactory.h src/individual.h
src/LCEcomposite.o: src/ttrait.h src/ttrait_with_map.h src/MPStatHandler.h
src/LCEcomposite.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/LCEcomposite.o: src/paramsparser.h src/MPImanager.h src/Uniform.h
src/LCEcomposite.o: src/LCEdisperse.h src/LCEselection.h
src/LCEdisperse.o: src/LCEdisperse.h src/lifecycleevent.h src/param.h
src/LCEdisperse.o: src/handler.h src/tmatrix.h src/output.h src/types.h
src/LCEdisperse.o: src/simcomponent.h src/fileservices.h src/service.h
src/LCEdisperse.o: src/statservices.h src/statrecorder.h
src/LCEdisperse.o: src/updaterservices.h src/binarystoragebuffer.h
src/LCEdisperse.o: src/metapop.h src/indfactory.h src/individual.h
src/LCEdisperse.o: src/ttrait.h src/ttrait_with_map.h src/MPStatHandler.h
src/LCEdisperse.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/LCEdisperse.o: src/paramsparser.h src/MPImanager.h src/Uniform.h
src/LCEmisc.o: src/LCEmisc.h src/types.h src/lifecycleevent.h src/param.h
src/LCEmisc.o: src/handler.h src/tmatrix.h src/output.h src/simcomponent.h
src/LCEmisc.o: src/fileservices.h src/service.h src/statservices.h
src/LCEmisc.o: src/statrecorder.h src/updaterservices.h
src/LCEmisc.o: src/binarystoragebuffer.h src/metapop.h src/indfactory.h
src/LCEmisc.o: src/individual.h src/ttrait.h src/ttrait_with_map.h
src/LCEmisc.o: src/MPStatHandler.h src/stathandler.h src/binarydataloader.h
src/LCEmisc.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/LCEmisc.o: src/filehandler.h src/Uniform.h src/tstring.h
src/LCEquanti.o: src/LCEquanti.h src/lifecycleevent.h src/param.h
src/LCEquanti.o: src/handler.h src/tmatrix.h src/output.h src/types.h
src/LCEquanti.o: src/simcomponent.h src/fileservices.h src/service.h
src/LCEquanti.o: src/statservices.h src/statrecorder.h src/updaterservices.h
src/LCEquanti.o: src/binarystoragebuffer.h src/metapop.h src/indfactory.h
src/LCEquanti.o: src/individual.h src/ttrait.h src/ttrait_with_map.h
src/LCEquanti.o: src/MPStatHandler.h src/stathandler.h src/binarydataloader.h
src/LCEquanti.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/LCEquanti.o: src/Uniform.h src/ttquanti.h src/filehandler.h
src/LCEquanti.o: src/datatable.h src/ttneutralgenes.h src/utils.h
src/LCEselection.o: src/LCEselection.h src/lifecycleevent.h src/param.h
src/LCEselection.o: src/handler.h src/tmatrix.h src/output.h src/types.h
src/LCEselection.o: src/simcomponent.h src/fileservices.h src/service.h
src/LCEselection.o: src/statservices.h src/statrecorder.h
src/LCEselection.o: src/updaterservices.h src/binarystoragebuffer.h
src/LCEselection.o: src/metapop.h src/indfactory.h src/individual.h
src/LCEselection.o: src/ttrait.h src/ttrait_with_map.h src/MPStatHandler.h
src/LCEselection.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/LCEselection.o: src/paramsparser.h src/MPImanager.h src/Uniform.h
src/LCEselection.o: src/utils.h src/tstring.h src/simenv.h src/simulation.h
src/LCEselection.o: src/basicsimulation.h
src/MPImanager.o: src/simulation.h src/basicsimulation.h src/ttrait.h
src/MPImanager.o: src/types.h src/simcomponent.h src/fileservices.h
src/MPImanager.o: src/service.h src/handler.h src/param.h src/tmatrix.h
src/MPImanager.o: src/output.h src/statservices.h src/statrecorder.h
src/MPImanager.o: src/updaterservices.h src/binarystoragebuffer.h
src/MPImanager.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/MPImanager.o: src/individual.h src/ttrait_with_map.h src/MPStatHandler.h
src/MPImanager.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/MPImanager.o: src/paramsparser.h src/MPImanager.h src/binarydatasaver.h
src/MPImanager.o: src/filehandler.h
src/MPStatHandler.o: src/MPStatHandler.h src/stathandler.h src/handler.h
src/MPStatHandler.o: src/statservices.h src/service.h src/statrecorder.h
src/MPStatHandler.o: src/types.h src/output.h src/metapop.h src/indfactory.h
src/MPStatHandler.o: src/individual.h src/ttrait.h src/simcomponent.h
src/MPStatHandler.o: src/fileservices.h src/param.h src/tmatrix.h
src/MPStatHandler.o: src/updaterservices.h src/binarystoragebuffer.h
src/MPStatHandler.o: src/ttrait_with_map.h src/binarydataloader.h
src/MPStatHandler.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/MPStatHandler.o: src/tstring.h
src/basicsimulation.o: src/basicsimulation.h src/ttrait.h src/types.h
src/basicsimulation.o: src/simcomponent.h src/fileservices.h src/service.h
src/basicsimulation.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/basicsimulation.o: src/statservices.h src/statrecorder.h
src/basicsimulation.o: src/updaterservices.h src/binarystoragebuffer.h
src/basicsimulation.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/basicsimulation.o: src/individual.h src/ttrait_with_map.h
src/basicsimulation.o: src/MPStatHandler.h src/stathandler.h
src/basicsimulation.o: src/binarydataloader.h src/fileparser.h
src/basicsimulation.o: src/paramsparser.h src/MPImanager.h
src/binarydataloader.o: src/tstring.h src/output.h src/Uniform.h
src/binarydataloader.o: src/binarydataloader.h src/binarystoragebuffer.h
src/binarydataloader.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/binarydataloader.o: src/basicsimulation.h src/ttrait.h src/types.h
src/binarydataloader.o: src/simcomponent.h src/fileservices.h src/service.h
src/binarydataloader.o: src/handler.h src/param.h src/tmatrix.h
src/binarydataloader.o: src/statservices.h src/statrecorder.h
src/binarydataloader.o: src/updaterservices.h src/lifecycleevent.h
src/binarydataloader.o: src/metapop.h src/indfactory.h src/individual.h
src/binarydataloader.o: src/ttrait_with_map.h src/MPStatHandler.h
src/binarydataloader.o: src/stathandler.h
src/binarydatasaver.o: src/metapop.h src/types.h src/indfactory.h
src/binarydatasaver.o: src/individual.h src/ttrait.h src/simcomponent.h
src/binarydatasaver.o: src/fileservices.h src/service.h src/handler.h
src/binarydatasaver.o: src/param.h src/tmatrix.h src/output.h
src/binarydatasaver.o: src/statservices.h src/statrecorder.h
src/binarydatasaver.o: src/updaterservices.h src/binarystoragebuffer.h
src/binarydatasaver.o: src/ttrait_with_map.h src/MPStatHandler.h
src/binarydatasaver.o: src/stathandler.h src/binarydataloader.h
src/binarydatasaver.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/binarydatasaver.o: src/binarydatasaver.h src/lifecycleevent.h
src/binarydatasaver.o: src/filehandler.h src/version.h src/simenv.h
src/binarydatasaver.o: src/simulation.h src/basicsimulation.h
src/bitstring.o: src/bitstring.h
src/filehandler.o: src/output.h src/filehandler.h src/handler.h
src/filehandler.o: src/fileservices.h src/service.h src/param.h src/tmatrix.h
src/filehandler.o: src/types.h src/LCEmisc.h src/lifecycleevent.h
src/filehandler.o: src/simcomponent.h src/statservices.h src/statrecorder.h
src/filehandler.o: src/updaterservices.h src/binarystoragebuffer.h
src/filehandler.o: src/metapop.h src/indfactory.h src/individual.h
src/filehandler.o: src/ttrait.h src/ttrait_with_map.h src/MPStatHandler.h
src/filehandler.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/filehandler.o: src/paramsparser.h src/MPImanager.h src/Uniform.h
src/filehandler.o: src/version.h
src/fileparser.o: src/fileparser.h src/paramsparser.h src/output.h
src/fileservices.o: src/simcomponent.h src/fileservices.h src/service.h
src/fileservices.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/fileservices.o: src/types.h src/statservices.h src/statrecorder.h
src/fileservices.o: src/updaterservices.h src/binarystoragebuffer.h
src/fileservices.o: src/filehandler.h src/metapop.h src/indfactory.h
src/fileservices.o: src/individual.h src/ttrait.h src/ttrait_with_map.h
src/fileservices.o: src/MPStatHandler.h src/stathandler.h
src/fileservices.o: src/binarydataloader.h src/fileparser.h
src/fileservices.o: src/paramsparser.h src/MPImanager.h
src/indfactory.o: src/indfactory.h src/individual.h src/types.h src/ttrait.h
src/indfactory.o: src/simcomponent.h src/fileservices.h src/service.h
src/indfactory.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/indfactory.o: src/statservices.h src/statrecorder.h src/updaterservices.h
src/indfactory.o: src/binarystoragebuffer.h src/ttrait_with_map.h
src/indfactory.o: src/Uniform.h
src/individual.o: src/individual.h src/types.h src/ttrait.h
src/individual.o: src/simcomponent.h src/fileservices.h src/service.h
src/individual.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/individual.o: src/statservices.h src/statrecorder.h src/updaterservices.h
src/individual.o: src/binarystoragebuffer.h src/ttrait_with_map.h
src/main.o: src/simulation.h src/basicsimulation.h src/ttrait.h src/types.h
src/main.o: src/simcomponent.h src/fileservices.h src/service.h src/handler.h
src/main.o: src/param.h src/tmatrix.h src/output.h src/statservices.h
src/main.o: src/statrecorder.h src/updaterservices.h
src/main.o: src/binarystoragebuffer.h src/lifecycleevent.h src/metapop.h
src/main.o: src/indfactory.h src/individual.h src/ttrait_with_map.h
src/main.o: src/MPStatHandler.h src/stathandler.h src/binarydataloader.h
src/main.o: src/fileparser.h src/paramsparser.h src/MPImanager.h src/simenv.h
src/metapop.o: src/Uniform.h src/output.h src/metapop.h src/types.h
src/metapop.o: src/indfactory.h src/individual.h src/ttrait.h
src/metapop.o: src/simcomponent.h src/fileservices.h src/service.h
src/metapop.o: src/handler.h src/param.h src/tmatrix.h src/statservices.h
src/metapop.o: src/statrecorder.h src/updaterservices.h
src/metapop.o: src/binarystoragebuffer.h src/ttrait_with_map.h
src/metapop.o: src/MPStatHandler.h src/stathandler.h src/binarydataloader.h
src/metapop.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/metapop.o: src/lifecycleevent.h src/simulation.h src/basicsimulation.h
src/metapop.o: src/simenv.h src/filehandler.h
src/output.o: src/output.h src/MPImanager.h
src/param.o: src/paramsparser.h src/param.h src/handler.h src/tmatrix.h
src/param.o: src/output.h src/types.h src/tstring.h
src/paramsparser.o: src/paramsparser.h src/output.h
src/patch.o: src/metapop.h src/types.h src/indfactory.h src/individual.h
src/patch.o: src/ttrait.h src/simcomponent.h src/fileservices.h src/service.h
src/patch.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/patch.o: src/statservices.h src/statrecorder.h src/updaterservices.h
src/patch.o: src/binarystoragebuffer.h src/ttrait_with_map.h
src/patch.o: src/MPStatHandler.h src/stathandler.h src/binarydataloader.h
src/patch.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/patch.o: src/Uniform.h
src/servicenotifiers.o: src/servicenotifiers.h src/types.h
src/servicenotifiers.o: src/lifecycleevent.h src/param.h src/handler.h
src/servicenotifiers.o: src/tmatrix.h src/output.h src/simcomponent.h
src/servicenotifiers.o: src/fileservices.h src/service.h src/statservices.h
src/servicenotifiers.o: src/statrecorder.h src/updaterservices.h
src/servicenotifiers.o: src/binarystoragebuffer.h src/metapop.h
src/servicenotifiers.o: src/indfactory.h src/individual.h src/ttrait.h
src/servicenotifiers.o: src/ttrait_with_map.h src/MPStatHandler.h
src/servicenotifiers.o: src/stathandler.h src/binarydataloader.h
src/servicenotifiers.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/servicenotifiers.o: src/filehandler.h src/Uniform.h src/simenv.h
src/servicenotifiers.o: src/simulation.h src/basicsimulation.h src/tstring.h
src/simenv.o: src/simenv.h src/simulation.h src/basicsimulation.h
src/simenv.o: src/ttrait.h src/types.h src/simcomponent.h src/fileservices.h
src/simenv.o: src/service.h src/handler.h src/param.h src/tmatrix.h
src/simenv.o: src/output.h src/statservices.h src/statrecorder.h
src/simenv.o: src/updaterservices.h src/binarystoragebuffer.h
src/simenv.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/simenv.o: src/individual.h src/ttrait_with_map.h src/MPStatHandler.h
src/simenv.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/simenv.o: src/paramsparser.h src/MPImanager.h src/binarydatasaver.h
src/simenv.o: src/filehandler.h src/LCEmisc.h src/Uniform.h src/LCEbreed.h
src/simenv.o: src/LCEdisperse.h src/LCEselection.h src/LCEcomposite.h
src/simenv.o: src/LCEquanti.h src/servicenotifiers.h src/ttneutralgenes.h
src/simenv.o: src/datatable.h src/ttdeletmutations_bitstring.h
src/simenv.o: src/bitstring.h src/ttdispersal.h src/ttwolbachia.h
src/simenv.o: src/ttquanti.h src/ttbdmi.h
src/simulation.o: src/simulation.h src/basicsimulation.h src/ttrait.h
src/simulation.o: src/types.h src/simcomponent.h src/fileservices.h
src/simulation.o: src/service.h src/handler.h src/param.h src/tmatrix.h
src/simulation.o: src/output.h src/statservices.h src/statrecorder.h
src/simulation.o: src/updaterservices.h src/binarystoragebuffer.h
src/simulation.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/simulation.o: src/individual.h src/ttrait_with_map.h src/MPStatHandler.h
src/simulation.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/simulation.o: src/paramsparser.h src/MPImanager.h src/servicenotifiers.h
src/simulation.o: src/filehandler.h src/Uniform.h src/version.h src/tstring.h
src/stathandler.o: src/stathandler.h src/handler.h src/statservices.h
src/stathandler.o: src/service.h src/statrecorder.h src/types.h src/output.h
src/stathandler.o: src/metapop.h src/indfactory.h src/individual.h
src/stathandler.o: src/ttrait.h src/simcomponent.h src/fileservices.h
src/stathandler.o: src/param.h src/tmatrix.h src/updaterservices.h
src/stathandler.o: src/binarystoragebuffer.h src/ttrait_with_map.h
src/stathandler.o: src/MPStatHandler.h src/binarydataloader.h
src/stathandler.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/stathandler.o: src/simenv.h src/simulation.h src/basicsimulation.h
src/stathandler.o: src/lifecycleevent.h
src/stats_coa.o: src/ttneutralgenes.h src/ttrait_with_map.h src/ttrait.h
src/stats_coa.o: src/types.h src/simcomponent.h src/fileservices.h
src/stats_coa.o: src/service.h src/handler.h src/param.h src/tmatrix.h
src/stats_coa.o: src/output.h src/statservices.h src/statrecorder.h
src/stats_coa.o: src/updaterservices.h src/binarystoragebuffer.h
src/stats_coa.o: src/filehandler.h src/stathandler.h src/datatable.h
src/stats_coa.o: src/metapop.h src/indfactory.h src/individual.h
src/stats_coa.o: src/MPStatHandler.h src/binarydataloader.h src/fileparser.h
src/stats_coa.o: src/paramsparser.h src/MPImanager.h
src/stats_delet_bitstring.o: src/ttdeletmutations_bitstring.h src/ttrait.h
src/stats_delet_bitstring.o: src/types.h src/simcomponent.h
src/stats_delet_bitstring.o: src/fileservices.h src/service.h src/handler.h
src/stats_delet_bitstring.o: src/param.h src/tmatrix.h src/output.h
src/stats_delet_bitstring.o: src/statservices.h src/statrecorder.h
src/stats_delet_bitstring.o: src/updaterservices.h src/binarystoragebuffer.h
src/stats_delet_bitstring.o: src/stathandler.h src/filehandler.h
src/stats_delet_bitstring.o: src/datatable.h src/metapop.h src/indfactory.h
src/stats_delet_bitstring.o: src/individual.h src/ttrait_with_map.h
src/stats_delet_bitstring.o: src/MPStatHandler.h src/binarydataloader.h
src/stats_delet_bitstring.o: src/fileparser.h src/paramsparser.h
src/stats_delet_bitstring.o: src/MPImanager.h src/bitstring.h src/Uniform.h
src/stats_demo.o: src/MPStatHandler.h src/stathandler.h src/handler.h
src/stats_demo.o: src/statservices.h src/service.h src/statrecorder.h
src/stats_demo.o: src/types.h src/output.h src/metapop.h src/indfactory.h
src/stats_demo.o: src/individual.h src/ttrait.h src/simcomponent.h
src/stats_demo.o: src/fileservices.h src/param.h src/tmatrix.h
src/stats_demo.o: src/updaterservices.h src/binarystoragebuffer.h
src/stats_demo.o: src/ttrait_with_map.h src/binarydataloader.h
src/stats_demo.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/stats_disp.o: src/ttdispersal.h src/ttrait.h src/types.h
src/stats_disp.o: src/simcomponent.h src/fileservices.h src/service.h
src/stats_disp.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/stats_disp.o: src/statservices.h src/statrecorder.h src/updaterservices.h
src/stats_disp.o: src/binarystoragebuffer.h src/stathandler.h src/metapop.h
src/stats_disp.o: src/indfactory.h src/individual.h src/ttrait_with_map.h
src/stats_disp.o: src/MPStatHandler.h src/binarydataloader.h src/fileparser.h
src/stats_disp.o: src/paramsparser.h src/MPImanager.h
src/stats_fstat.o: src/metapop.h src/types.h src/indfactory.h
src/stats_fstat.o: src/individual.h src/ttrait.h src/simcomponent.h
src/stats_fstat.o: src/fileservices.h src/service.h src/handler.h src/param.h
src/stats_fstat.o: src/tmatrix.h src/output.h src/statservices.h
src/stats_fstat.o: src/statrecorder.h src/updaterservices.h
src/stats_fstat.o: src/binarystoragebuffer.h src/ttrait_with_map.h
src/stats_fstat.o: src/MPStatHandler.h src/stathandler.h
src/stats_fstat.o: src/binarydataloader.h src/fileparser.h src/paramsparser.h
src/stats_fstat.o: src/MPImanager.h src/ttneutralgenes.h src/filehandler.h
src/stats_fstat.o: src/datatable.h
src/stats_wolbachia.o: src/ttwolbachia.h src/ttrait.h src/types.h
src/stats_wolbachia.o: src/simcomponent.h src/fileservices.h src/service.h
src/stats_wolbachia.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/stats_wolbachia.o: src/statservices.h src/statrecorder.h
src/stats_wolbachia.o: src/updaterservices.h src/binarystoragebuffer.h
src/stats_wolbachia.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/stats_wolbachia.o: src/individual.h src/ttrait_with_map.h
src/stats_wolbachia.o: src/MPStatHandler.h src/stathandler.h
src/stats_wolbachia.o: src/binarydataloader.h src/fileparser.h
src/stats_wolbachia.o: src/paramsparser.h src/MPImanager.h src/filehandler.h
src/stats_wolbachia.o: src/Uniform.h src/LCEbreed.h
src/statservices.o: src/statservices.h src/service.h src/handler.h
src/statservices.o: src/statrecorder.h src/types.h src/output.h
src/statservices.o: src/stathandler.h src/metapop.h src/indfactory.h
src/statservices.o: src/individual.h src/ttrait.h src/simcomponent.h
src/statservices.o: src/fileservices.h src/param.h src/tmatrix.h
src/statservices.o: src/updaterservices.h src/binarystoragebuffer.h
src/statservices.o: src/ttrait_with_map.h src/MPStatHandler.h
src/statservices.o: src/binarydataloader.h src/fileparser.h
src/statservices.o: src/paramsparser.h src/MPImanager.h src/simenv.h
src/statservices.o: src/simulation.h src/basicsimulation.h
src/statservices.o: src/lifecycleevent.h
src/ttbdmi.o: src/ttbdmi.h src/ttrait.h src/types.h src/simcomponent.h
src/ttbdmi.o: src/fileservices.h src/service.h src/handler.h src/param.h
src/ttbdmi.o: src/tmatrix.h src/output.h src/statservices.h
src/ttbdmi.o: src/statrecorder.h src/updaterservices.h
src/ttbdmi.o: src/binarystoragebuffer.h src/stathandler.h src/filehandler.h
src/ttbdmi.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/ttbdmi.o: src/individual.h src/ttrait_with_map.h src/MPStatHandler.h
src/ttbdmi.o: src/binarydataloader.h src/fileparser.h src/paramsparser.h
src/ttbdmi.o: src/MPImanager.h src/bitstring.h src/Uniform.h
src/ttdeletmutations_bitstring.o: src/ttdeletmutations_bitstring.h
src/ttdeletmutations_bitstring.o: src/ttrait.h src/types.h src/simcomponent.h
src/ttdeletmutations_bitstring.o: src/fileservices.h src/service.h
src/ttdeletmutations_bitstring.o: src/handler.h src/param.h src/tmatrix.h
src/ttdeletmutations_bitstring.o: src/output.h src/statservices.h
src/ttdeletmutations_bitstring.o: src/statrecorder.h src/updaterservices.h
src/ttdeletmutations_bitstring.o: src/binarystoragebuffer.h src/stathandler.h
src/ttdeletmutations_bitstring.o: src/filehandler.h src/datatable.h
src/ttdeletmutations_bitstring.o: src/metapop.h src/indfactory.h
src/ttdeletmutations_bitstring.o: src/individual.h src/ttrait_with_map.h
src/ttdeletmutations_bitstring.o: src/MPStatHandler.h src/binarydataloader.h
src/ttdeletmutations_bitstring.o: src/fileparser.h src/paramsparser.h
src/ttdeletmutations_bitstring.o: src/MPImanager.h src/bitstring.h
src/ttdeletmutations_bitstring.o: src/Uniform.h
src/ttdispersal.o: src/ttdispersal.h src/ttrait.h src/types.h
src/ttdispersal.o: src/simcomponent.h src/fileservices.h src/service.h
src/ttdispersal.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/ttdispersal.o: src/statservices.h src/statrecorder.h
src/ttdispersal.o: src/updaterservices.h src/binarystoragebuffer.h
src/ttdispersal.o: src/stathandler.h src/metapop.h src/indfactory.h
src/ttdispersal.o: src/individual.h src/ttrait_with_map.h src/MPStatHandler.h
src/ttdispersal.o: src/binarydataloader.h src/fileparser.h src/paramsparser.h
src/ttdispersal.o: src/MPImanager.h src/Uniform.h
src/ttneutralgenes.o: src/ttneutralgenes.h src/ttrait_with_map.h src/ttrait.h
src/ttneutralgenes.o: src/types.h src/simcomponent.h src/fileservices.h
src/ttneutralgenes.o: src/service.h src/handler.h src/param.h src/tmatrix.h
src/ttneutralgenes.o: src/output.h src/statservices.h src/statrecorder.h
src/ttneutralgenes.o: src/updaterservices.h src/binarystoragebuffer.h
src/ttneutralgenes.o: src/filehandler.h src/stathandler.h src/datatable.h
src/ttneutralgenes.o: src/metapop.h src/indfactory.h src/individual.h
src/ttneutralgenes.o: src/MPStatHandler.h src/binarydataloader.h
src/ttneutralgenes.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/ttneutralgenes.o: src/Uniform.h src/tstring.h
src/ttquanti.o: src/ttquanti.h src/ttrait_with_map.h src/ttrait.h src/types.h
src/ttquanti.o: src/simcomponent.h src/fileservices.h src/service.h
src/ttquanti.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/ttquanti.o: src/statservices.h src/statrecorder.h src/updaterservices.h
src/ttquanti.o: src/binarystoragebuffer.h src/filehandler.h src/stathandler.h
src/ttquanti.o: src/metapop.h src/indfactory.h src/individual.h
src/ttquanti.o: src/MPStatHandler.h src/binarydataloader.h src/fileparser.h
src/ttquanti.o: src/paramsparser.h src/MPImanager.h src/datatable.h
src/ttquanti.o: src/Uniform.h src/tstring.h src/utils.h
src/ttrait_with_map.o: src/ttrait_with_map.h src/ttrait.h src/types.h
src/ttrait_with_map.o: src/simcomponent.h src/fileservices.h src/service.h
src/ttrait_with_map.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/ttrait_with_map.o: src/statservices.h src/statrecorder.h
src/ttrait_with_map.o: src/updaterservices.h src/binarystoragebuffer.h
src/ttrait_with_map.o: src/Uniform.h src/tstring.h
src/ttwolbachia.o: src/ttwolbachia.h src/ttrait.h src/types.h
src/ttwolbachia.o: src/simcomponent.h src/fileservices.h src/service.h
src/ttwolbachia.o: src/handler.h src/param.h src/tmatrix.h src/output.h
src/ttwolbachia.o: src/statservices.h src/statrecorder.h
src/ttwolbachia.o: src/updaterservices.h src/binarystoragebuffer.h
src/ttwolbachia.o: src/lifecycleevent.h src/metapop.h src/indfactory.h
src/ttwolbachia.o: src/individual.h src/ttrait_with_map.h src/MPStatHandler.h
src/ttwolbachia.o: src/stathandler.h src/binarydataloader.h src/fileparser.h
src/ttwolbachia.o: src/paramsparser.h src/MPImanager.h src/filehandler.h
src/ttwolbachia.o: src/Uniform.h src/LCEbreed.h src/simenv.h src/simulation.h
src/ttwolbachia.o: src/basicsimulation.h
src/updaterservices.o: src/updaterservices.h src/service.h src/handler.h
src/updaterservices.o: src/param.h src/tmatrix.h src/output.h src/types.h
src/updaterservices.o: src/simcomponent.h src/fileservices.h
src/updaterservices.o: src/statservices.h src/statrecorder.h
src/updaterservices.o: src/binarystoragebuffer.h src/lifecycleevent.h
src/updaterservices.o: src/metapop.h src/indfactory.h src/individual.h
src/updaterservices.o: src/ttrait.h src/ttrait_with_map.h src/MPStatHandler.h
src/updaterservices.o: src/stathandler.h src/binarydataloader.h
src/updaterservices.o: src/fileparser.h src/paramsparser.h src/MPImanager.h
src/utils.o: src/utils.h src/types.h src/tmatrix.h src/output.h src/Uniform.h
