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
	${OBJECTDIR}/2d_malloc.o \
	${OBJECTDIR}/file_io.o \
	${OBJECTDIR}/filter.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/main_async.o \
	${OBJECTDIR}/main_async_nonper.o \
	${OBJECTDIR}/main_async_omp.o \
	${OBJECTDIR}/main_async_omp_simple.o \
	${OBJECTDIR}/main_serial.o \
	${OBJECTDIR}/main_serial_omp.o \
	${OBJECTDIR}/main_sync.o \
	${OBJECTDIR}/main_sync_omp.o \
	${OBJECTDIR}/main_sync_omp_simple.o \
	${OBJECTDIR}/topology.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/par_conv

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/par_conv: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.c} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/par_conv ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/2d_malloc.o: 2d_malloc.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/2d_malloc.o 2d_malloc.c

${OBJECTDIR}/file_io.o: file_io.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/file_io.o file_io.c

${OBJECTDIR}/filter.o: filter.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/filter.o filter.c

${OBJECTDIR}/main.o: main.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.c

${OBJECTDIR}/main_async.o: main_async.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_async.o main_async.c

${OBJECTDIR}/main_async_nonper.o: main_async_nonper.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_async_nonper.o main_async_nonper.c

${OBJECTDIR}/main_async_omp.o: main_async_omp.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_async_omp.o main_async_omp.c

${OBJECTDIR}/main_async_omp_simple.o: main_async_omp_simple.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_async_omp_simple.o main_async_omp_simple.c

${OBJECTDIR}/main_serial.o: main_serial.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_serial.o main_serial.c

${OBJECTDIR}/main_serial_omp.o: main_serial_omp.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_serial_omp.o main_serial_omp.c

${OBJECTDIR}/main_sync.o: main_sync.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_sync.o main_sync.c

${OBJECTDIR}/main_sync_omp.o: main_sync_omp.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_sync_omp.o main_sync_omp.c

${OBJECTDIR}/main_sync_omp_simple.o: main_sync_omp_simple.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_sync_omp_simple.o main_sync_omp_simple.c

${OBJECTDIR}/topology.o: topology.c 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -Wall -I/usr/local/mpich2/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/topology.o topology.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/par_conv

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
