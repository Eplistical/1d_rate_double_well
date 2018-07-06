# - Try to find the Intel Math Kernel Library
# Once done this will define
#
#  MKL_FOUND - system has MKL
#  MKL_ROOT_DIR - path to the MKL base directory
#  MKL_INCLUDE_DIRS - the MKL INCLUDE directory
#  MKL_LIBRARIES - MKL libraries
#
# There are few SETs of libraries:

# Array indexes modes:
#   LP - 32 bit indexes of arrays
#   ILP - 64 bit indexes of arrays
# Threading:
#   SEQUENTIAL - no threading
#   INTEL - Intel threading library
#   GNU - GNU threading library
# MPI support
#   NOMPI - no MPI support
#   INTEL - Intel MPI library
#   OPEN - Open MPI library
#   SGI - SGI MPT Library

# architecture
IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    SET(MKL_ARCH 64)
    SET(MKL_ARCH_DIR "intel64")
ELSE()
    SET(MKL_ARCH 32)
    SET(MKL_ARCH_DIR "ia32")
ENDIF()

SET(MKL_THREAD_VARIANTS SEQUENTIAL GNUTHREAD INTELTHREAD)
SET(MKL_MODE_VARIANTS ILP LP)
SET(MKL_MPI_VARIANTS NOMPI INTELMPI OPENMPI SGIMPT)

SET(CMAKE_FIND_DEBUG_MODE 1)

SET(MKL_ROOT_DIR_CANDIDATES
    $ENV{MKLROOT}
    $ENV{MKL_ROOT}
    $ENV{MKLDIR}
    $ENV{MKL_DIR}
    /opt/intel/mkl
    /opt/intel/cmkl
    /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
    "C:/Program Files (x86)/Intel/ComposerXE-2011/mkl"
    "C:/Program Files (x86)/Intel/Composer XE 2013/mkl"
    "C:/Program Files/Intel/MKL/*/"
    "C:/Program Files/Intel/ComposerXE-2011/mkl"
    "C:/Program Files/Intel/Composer XE 2013/mkl"
    "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/"
    "C:/Program Files/Intel/Composer XE 2015/mkl/"
    )

FIND_PATH(MKL_ROOT_DIR NAMES include/mkl_cblas.h PATHS ${MKL_ROOT_DIR_CANDIDATES})

IF(MKL_ROOT_DIR)
    MESSAGE("-- MKL found at ${MKL_ROOT_DIR}")
    SET(HAVE_MKL TRUE)
ELSE(MKL_ROOT_DIR)
    MESSAGE(WARNING "MKL not found")
    SET(HAVE_MKL FALSE)
ENDIF(MKL_ROOT_DIR)

IF(HAVE_MKL)
    FIND_PATH(MKL_INCLUDE_DIRS mkl_cblas.h PATHS ${MKL_ROOT_DIR}/include)

    FIND_PATH(MKL_FFTW_INCLUDE_DIR fftw3.h PATH_SUFFIXES fftw PATHS ${MKL_ROOT_DIR}/include NO_DEFAULT_PATH)

    IF(WIN32)
        SET(MKL_LIB_SEARCHPATH $ENV{ICC_LIB_DIR} $ENV{MKL_LIB_DIR} "${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}" "${MKL_ROOT_DIR}/../compiler" "${MKL_ROOT_DIR}/../compiler/lib/${MKL_ARCH_DIR}")

        IF(MKL_INCLUDE_DIRS MATCHES "10.")
            SET(MKL_LIBS mkl_solver mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libguide mkl_lapack95 mkl_blas95)
            IF(${MKL_ARCH} EQUAL 64)
                SET(MKL_LIBS mkl_solver_lp64 mkl_core mkl_intel_lp64 mkl_intel_thread libguide mkl_lapack95_lp64 mkl_blas95_lp64)
            ENDIF()
        ELSEIF(MKL_INCLUDE_DIRS MATCHES "2013") # version 11 ...
            SET(MKL_LIBS mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libiomp5md mkl_lapack95 mkl_blas95)
            IF(${MKL_ARCH} EQUAL 64)
                SET(MKL_LIBS mkl_core mkl_intel_lp64 mkl_intel_thread libiomp5md mkl_lapack95_lp64 mkl_blas95_lp64)
            ENDIF()
        ELSEIF(MKL_INCLUDE_DIRS MATCHES "2015")
            IF(${MKL_ARCH} EQUAL 64)
                SET(MKL_LIBS mkl_intel_lp64 mkl_core mkl_intel_thread mkl_lapack95_lp64 mkl_blas95_lp64 )
            ELSE()
                SET(MKL_LIBS mkl_intel_c mkl_core mkl_intel_thread mkl_lapack95 mkl_blas95 )
            ENDIF()
        ELSE() # old MKL 9
            SET(MKL_LIBS mkl_solver mkl_c libguide mkl_lapack mkl_ia32)
        ENDIF()

        IF(MKL_INCLUDE_DIRS MATCHES "10.3")
            SET(MKL_LIBS ${MKL_LIBS} libiomp5md)
        ENDIF()

        FOREACH (LIB ${MKL_LIBS})
            FIND_LIBRARY(${LIB}_PATH ${LIB} PATHS ${MKL_LIB_SEARCHPATH} ENV LIBRARY_PATH)
            IF(${LIB}_PATH)
                SET(MKL_LIBRARIES ${MKL_LIBRARIES} ${${LIB}_PATH})
            ELSE()
                MESSAGE(FATAL_ERROR "MKL: cannot find ${LIB}")
                BREAK()
            ENDIF()
        ENDFOREACH()
        SET(MKL_FOUND ON)

    ELSE() # Unix-like system

        SET(MKL_LIB_DIR_CANDIDATES ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR} ${MKL_ROOT_DIR}/lib)

        FIND_LIBRARY(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_LIB_DIR_CANDIDATES})

        # Threading libraries

        FIND_LIBRARY(MKL_RT_LIBRARY mkl_rt PATHS ${MKL_LIB_DIR_CANDIDATES})
        FIND_LIBRARY(MKL_SEQUENTIAL_LIBRARY mkl_sequential PATHS ${MKL_LIB_DIR_CANDIDATES})
        FIND_LIBRARY(MKL_INTELTHREAD_LIBRARY mkl_intel_thread PATHS ${MKL_LIB_DIR_CANDIDATES})
        FIND_LIBRARY(MKL_GNUTHREAD_LIBRARY mkl_gnu_thread PATHS ${MKL_LIB_DIR_CANDIDATES})

        # Intel Libraries

        IF(${MKL_ARCH} EQUAL 64)
            SET(INTEL_LP_SUFFIX  "_lp64")
            SET(INTEL_ILP_SUFFIX "_ilp64")
        ENDIF()

        FIND_LIBRARY(MKL_LP_LIBRARY mkl_intel${INTEL_LP_SUFFIX} PATHS ${MKL_LIB_DIR_CANDIDATES})
        FIND_LIBRARY(MKL_ILP_LIBRARY mkl_intel${INTEL_ILP_SUFFIX} PATHS ${MKL_LIB_DIR_CANDIDATES})

        # Lapack

        FIND_LIBRARY(MKL_LAPACK_LIBRARY mkl_lapack PATHS ${MKL_LIB_DIR_CANDIDATES})

        IF(NOT MKL_LAPACK_LIBRARY)
            FIND_LIBRARY(MKL_LAPACK_LIBRARY mkl_lapack95_lp64 PATHS ${MKL_LIB_DIR_CANDIDATES})
        ENDIF()

        # iomp5

        IF(UNIX AND NOT APPLE)
            FIND_LIBRARY(MKL_IOMP5_LIBRARY iomp5 PATHS ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR})
        ENDIF()

        FOREACH (MODEVAR ${MKL_MODE_VARIANTS})
            FOREACH (THREADVAR ${MKL_THREAD_VARIANTS})
                IF(MKL_CORE_LIBRARY AND MKL_${MODEVAR}_LIBRARY AND MKL_${THREADVAR}_LIBRARY)
                    SET(MKL_${MODEVAR}_${THREADVAR}_LIBRARIES
                        ${MKL_${MODEVAR}_LIBRARY} ${MKL_${THREADVAR}_LIBRARY} ${MKL_CORE_LIBRARY}
                        ${MKL_LAPACK_LIBRARY} ${MKL_IOMP5_LIBRARY})
                    #MESSAGE("${MODEVAR} ${THREADVAR} ${MKL_${MODEVAR}_${THREADVAR}_LIBRARIES}") # for debug
                ENDIF()
            ENDFOREACH()
        ENDFOREACH()

        SET(MKL_LIBRARIES ${MKL_RT_LIBRARY})
        MARK_AS_ADVANCED(MKL_CORE_LIBRARY MKL_LP_LIBRARY MKL_ILP_LIBRARY
            MKL_SEQUENTIAL_LIBRARY MKL_INTELTHREAD_LIBRARY MKL_GNUTHREAD_LIBRARY)
    ENDIF()

    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_INCLUDE_DIRS MKL_LIBRARIES)

    MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES)
ENDIF()
