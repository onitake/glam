cmake_minimum_required(VERSION 2.6)

project(GLAM)

set(GLAM_VERSION_MAJOR 0)
set(GLAM_VERSION_MINOR 3)
set(GLAM_VERSION_MICRO 0)

# Run cmake -D CMAKE_BUILD_TYPE=Debug .
# or  cmake -D CMAKE_BUILD_TYPE=Release .
# to enable debugging or optimizations from the command line.
# Debug builds have range checks enabled by default.
if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
	set(CMAKE_BUILD_TYPE Debug)
	#set(CMAKE_BUILD_TYPE Release)
endif()

include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(CheckIncludeFileCXX)

find_package(CxxTest)
if(CXXTEST_FOUND)
	set(GLAM_HAS_CXXTEST 1)
	include_directories(${CXXTEST_INCLUDE_DIR})
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(GLAM_DEBUG 1)
	check_cxx_compiler_flag("-O0 -g" GLAM_HAS_CXX_OPTION_DEBUG_GCC)
	if(GLAM_HAS_CXX_OPTION_DEBUG_GCC)
		set(GLAM_DEBUG_FLAGS "-O0 -g")
	else()
		check_cxx_compiler_flag("/Od /Zi /DEBUG" GLAM_HAS_CXX_OPTION_DEBUG_VC)
		if(GLAM_HAS_CXX_OPTION_DEBUG_VC)
			set(GLAM_DEBUG_FLAGS "/Od /Zi /DEBUG")
		endif()
	endif()
	if(GLAM_DEBUG_FLAGS)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLAM_DEBUG_FLAGS}")
	else()
		message(WARNING "Unknown compiler. Don't know how to set debug flags.")
	endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
	unset(GLAM_DEBUG)
	set(GLAM_RANGE_CHECKS)
	check_cxx_compiler_flag("-Wall -Ofast -ffast-math -ftree-vectorize -ftree-parallelize-loops=4 -g" GLAM_HAS_CXX_OPTION_OPTIMIZE_GCC)
	if(GLAM_HAS_CXX_OPTION_OPTIMIZE_GCC)
		set(GLAM_OPTIMIZE_FLAGS "-Wall -Ofast -ffast-math -ftree-vectorize -ftree-parallelize-loops=4 -g")
	else()
		check_cxx_compiler_flag("/O2 /Qpar /Zi /DEBUG" GLAM_HAS_CXX_OPTION_OPTIMIZE_VC)
		if(GLAM_HAS_CXX_OPTION_OPTIMIZE_VC)
			set(GLAM_OPTIMIZE_FLAGS "/O2 /Qpar /Zi /DEBUG")
		endif()
	endif()
	if(GLAM_OPTIMIZE_FLAGS)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLAM_OPTIMIZE_FLAGS}")
	else()
		message(WARNING "Unknown compiler. Don't know how to set optimize flags.")
	endif()
else()
	message(FATAL_ERROR "Unknown build type ${CMAKE_BUILD_TYPE}.")
endif()

check_cxx_compiler_flag("-std=c++11" GLAM_HAS_CXX_OPTION_CXX11_GCC47)
if(GLAM_HAS_CXX_OPTION_CXX11_GCC47)
	set(GLAM_CXX11_FLAGS "-std=c++11")
else()
	check_cxx_compiler_flag("-std=c++0x" GLAM_HAS_CXX_OPTION_CXX11_GCC43)
	if(GLAM_HAS_CXX_OPTION_CXX11_GCC43)
		set(GLAM_CXX11_FLAGS "-std=c++0x")
	endif()
endif()
if(DEFINED GLAM_CXX11_FLAGS)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLAM_CXX11_FLAGS}")
endif()

check_cxx_source_compiles("
int main() {
	long double d = 1.0L;
	return 0;
}
" GLAM_HAS_LONG_DOUBLE)
if(GLAM_HAS_LONG_DOUBLE)
else()
	message(WARNING "Your compiler doesn't support long double.")
endif()

check_cxx_source_compiles("
int f(int &&x) { return x; }
int main() {
	f(10);
	return 0;
}
" GLAM_HAS_XVALUE)
if(GLAM_HAS_XVALUE)
else()
	message(WARNING "Your compiler doesn't support C++11 xvalues.")
endif()

check_cxx_source_compiles("
#include <tuple>
struct Test { };
int main() {
	std::tuple<int, float, Test> t = std::make_tuple(1, 0.5f, Test());
	return 0;
}
" GLAM_HAS_TUPLE)
if(GLAM_HAS_TUPLE)
else()
	message(WARNING "Your compiler doesn't support C++11 tuples.")
endif()

check_cxx_source_compiles("
int main() {
	static_assert(true, \"true is true\");
	return 0;
}
" GLAM_HAS_STATIC_ASSERT)
if(GLAM_HAS_STATIC_ASSERT)
else()
	message(FATAL_ERROR "Your compiler doesn't support C++11 static assert.")
endif()

check_cxx_source_compiles("
#include <cmath>
int main() {
	std::round(0.5);
	return 0;
}
" GLAM_HAS_CMATH)
if(GLAM_HAS_CMATH)
else()
	message(FATAL_ERROR "Your compiler doesn't support C++11 <cmath>.")
endif()

check_cxx_source_compiles("
void test(int i) { }
template <typename... Args>
void test(int i, Args... args) {
	test(args...);
}
int main() {
	test(1, 2, 3, 4, 5);
	return 0;
}
" GLAM_HAS_VARIADIC_TEMPLATES)
if(GLAM_HAS_VARIADIC_TEMPLATES)
else()
	message(FATAL_ERROR "Your compiler doesn't support C++11 variadic templates.")
endif()

check_cxx_source_compiles("
struct TEST {
	static constexpr int CON = 1;
};
int main() {
	return TEST::CON - 1;
}
" GLAM_HAS_CONSTEXPR)
if(GLAM_HAS_CONSTEXPR)
else()
	message(FATAL_ERROR "Your compiler doesn't support C++11 constexpr.")
endif()

find_path(GMP_INCLUDE_DIRS gmp.h)
find_library(GMP_LIBRARIES NAMES gmp)
find_path(GMPXX_INCLUDE_DIRS gmpxx.h)
find_library(GMPXX_LIBRARIES NAMES gmpxx)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIBRARIES GMP_INCLUDE_DIRS)
find_package_handle_standard_args(GMPXX DEFAULT_MSG GMPXX_LIBRARIES GMPXX_INCLUDE_DIRS)

configure_file(
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/glam/config.h"
)

enable_testing()

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(glam)
add_subdirectory(test)
add_subdirectory(example)
