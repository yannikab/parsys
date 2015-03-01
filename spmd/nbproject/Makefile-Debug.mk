#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=/usr/local/mpich2/bin/mpicc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/async/main_async.o \
	${OBJECTDIR}/async/main_async_nonper.o \
	${OBJECTDIR}/async/main_async_omp.o \
	${OBJECTDIR}/async/main_async_omp_simple.o \
	${OBJECTDIR}/common/2d_malloc.o \
	${OBJECTDIR}/common/file_io.o \
	${OBJECTDIR}/common/filter.o \
	${OBJECTDIR}/common/topology.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/serial/main_serial.o \
	${OBJECTDIR}/serial/main_serial_omp.o \
	${OBJECTDIR}/sync/main_sync.o \
	${OBJECTDIR}/sync/main_sync_omp.o \
	${OBJECTDIR}/sync/main_sync_omp_simple.o


# C Compiler Flags
CFLAGS=-mpe=mpilog -fopenmp

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spmd

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spmd: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spmd ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/async/main_async.o: async/main_async.c 
	${MKDIR} -p ${OBJECTDIR}/async
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/async/main_async.o async/main_async.c

${OBJECTDIR}/async/main_async_nonper.o: async/main_async_nonper.c 
	${MKDIR} -p ${OBJECTDIR}/async
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/async/main_async_nonper.o async/main_async_nonper.c

${OBJECTDIR}/async/main_async_omp.o: async/main_async_omp.c 
	${MKDIR} -p ${OBJECTDIR}/async
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/async/main_async_omp.o async/main_async_omp.c

${OBJECTDIR}/async/main_async_omp_simple.o: async/main_async_omp_simple.c 
	${MKDIR} -p ${OBJECTDIR}/async
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/async/main_async_omp_simple.o async/main_async_omp_simple.c

${OBJECTDIR}/common/2d_malloc.o: common/2d_malloc.c 
	${MKDIR} -p ${OBJECTDIR}/common
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/common/2d_malloc.o common/2d_malloc.c

${OBJECTDIR}/common/file_io.o: common/file_io.c 
	${MKDIR} -p ${OBJECTDIR}/common
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/common/file_io.o common/file_io.c

${OBJECTDIR}/common/filter.o: common/filter.c 
	${MKDIR} -p ${OBJECTDIR}/common
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/common/filter.o common/filter.c

${OBJECTDIR}/common/topology.o: common/topology.c 
	${MKDIR} -p ${OBJECTDIR}/common
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/common/topology.o common/topology.c

${OBJECTDIR}/main.o: main.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.c

${OBJECTDIR}/serial/main_serial.o: serial/main_serial.c 
	${MKDIR} -p ${OBJECTDIR}/serial
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/serial/main_serial.o serial/main_serial.c

${OBJECTDIR}/serial/main_serial_omp.o: serial/main_serial_omp.c 
	${MKDIR} -p ${OBJECTDIR}/serial
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/serial/main_serial_omp.o serial/main_serial_omp.c

${OBJECTDIR}/sync/main_sync.o: sync/main_sync.c 
	${MKDIR} -p ${OBJECTDIR}/sync
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sync/main_sync.o sync/main_sync.c

${OBJECTDIR}/sync/main_sync_omp.o: sync/main_sync_omp.c 
	${MKDIR} -p ${OBJECTDIR}/sync
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sync/main_sync_omp.o sync/main_sync_omp.c

${OBJECTDIR}/sync/main_sync_omp_simple.o: sync/main_sync_omp_simple.c 
	${MKDIR} -p ${OBJECTDIR}/sync
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sync/main_sync_omp_simple.o sync/main_sync_omp_simple.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spmd

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
