cmake_minimum_required(VERSION 3.20)
project(MyProject CXX)

# ----- Options matching Makefile conditionals -----
option(DEBUG "Enable debug build" OFF)
option(PROFILE "Enable profiling" OFF)
option(VALGRIND "Enable valgrind checks" OFF)
option(FAST_ERFC "Enable FAST_ERFC" OFF)
option(LATTICE_DYNAMICS "Enable LATTICE_DYNAMICS" OFF)

set(CMAKE_CXX_STANDARD 11)

# ----- Base flags -----
set(BASE_CFLAGS "-DSFMT_MEXP=19937 -DHAVE_SSE2")
set(BASE_CFLAGS "${BASE_CFLAGS} -msse2")

set(BASE_FLAGS "-Wall")

# Optimization/debug/profile flags
if (DEBUG)
    set(BASE_FLAGS "${BASE_FLAGS} -O0 -g -DDEBUG")
    if (VALGRIND)
        set(BASE_CFLAGS "${BASE_CFLAGS} -DVALGRIND_CHECKS")
    endif()
else()
    if (PROFILE)
        set(BASE_FLAGS "${BASE_FLAGS} -O1 -pg")
    else()
        set(BASE_FLAGS "${BASE_FLAGS} -O3")
    endif()
endif()

if (FAST_ERFC)
    set(BASE_CFLAGS "${BASE_CFLAGS} -DFAST_ERFC")
endif()

if (LATTICE_DYNAMICS)
    set(BASE_CFLAGS "${BASE_CFLAGS} -DLATTICE_DYNAMICS")
    include_directories(/usr/include/eigen3)  # On Windows you'll change this path
endif()

# Apply flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BASE_FLAGS} ${BASE_CFLAGS}")

# ----- Source files -----
file(GLOB ALL_CPP *.cpp)
file(GLOB HEADERS *.h)

# Exclude main-* from general sources
set(SRC_FILES "")
foreach(f ${ALL_CPP})
    if (NOT f MATCHES "main-.*\\.cpp$")
        list(APPEND SRC_FILES ${f})
    endif()
endforeach()

# Mains
set(MAIN_NAMES virial-trappe)

# ----- Build each main-* executable -----
foreach(name ${MAIN_NAMES})
    set(main_src "main-${name}.cpp"
            TraPPEParams.cpp
            TraPPEParams.h
            unit.h
            constants.h
            P2PotentialGroupBuilder.cpp
            P2PotentialGroupBuilder.h
            unit.cpp
            move-cluster-angle-general.cpp
            move-cluster-stretch.cpp)
    add_executable(${name} ${SRC_FILES} ${main_src})

    # Windows doesn't need -lm; on Linux/macOS link against libm
    if (UNIX)
        target_link_libraries(${name} m)
    endif()
endforeach()

