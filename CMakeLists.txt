cmake_minimum_required(VERSION 2.6)

project(GLAM)

set(GLAM_VERSION_MAJOR 0)
set(GLAM_VERSION_MINOR 2)
set(GLAM_VERSION_MICRO 1)
set(GLAM_DEBUG 0)

include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(CheckIncludeFileCXX)

find_package(CxxTest)
if(CXXTEST_FOUND)
	set(HAS_CXXTEST 1)
	include_directories(${CXXTEST_INCLUDE_DIR})
endif()

if(GLAM_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
	check_cxx_compiler_flag("-O0 -g" HAS_CXX_OPTION_DEBUG_GCC)
	if (HAS_CXX_OPTION_DEBUG_GCC)
		set(GLAM_DEBUG_FLAGS "-O0 -g")
	else()
		check_cxx_compiler_flag("/Od /Zi /DEBUG" HAS_CXX_OPTION_DEBUG_VC)
		if(HAS_CXX_OPTION_DEBUG_VC)
			set(GLAM_DEBUG_FLAGS "/Od /Zi /DEBUG")
		endif()
	endif()
	if(GLAM_DEBUG_FLAGS)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLAM_DEBUG_FLAGS}")
	else()
		message(WARNING "Unknown compiler. Don't know how to set debug flags.")
	endif()
else()
	set(CMAKE_BUILD_TYPE Release)
	check_cxx_compiler_flag("-Wall -Ofast -ffast-math -ftree-vectorize -ftree-parallelize-loops=4 -g" HAS_CXX_OPTION_OPTIMIZE_GCC)
	if (HAS_CXX_OPTION_OPTIMIZE_GCC)
		set(GLAM_OPTIMIZE_FLAGS "-Wall -Ofast -ffast-math -ftree-vectorize -ftree-parallelize-loops=4 -g")
	else()
		check_cxx_compiler_flag("/O2 /Zi /DEBUG" HAS_CXX_OPTION_OPTIMIZE_VC)
		if(HAS_CXX_OPTION_OPTIMIZE_VC)
			set(GLAM_OPTIMIZE_FLAGS "/O2 /Zi /DEBUG")
		endif()
	endif()
	if(GLAM_OPTIMIZE_FLAGS)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLAM_OPTIMIZE_FLAGS}")
	else()
		message(WARNING "Unknown compiler. Don't know how to set optimize flags.")
	endif()
endif()

check_cxx_compiler_flag("-std=c++11" HAS_CXX_OPTION_CXX11_GCC47)
if (HAS_CXX_OPTION_CXX11_GCC47)
	set(GLAM_CXX11_FLAGS "-std=c++11")
else()
	check_cxx_compiler_flag("-std=c++0x" HAS_CXX_OPTION_CXX11_GCC43)
	if(HAS_CXX_OPTION_CXX11_GCC43)
		set(GLAM_CXX11_FLAGS "-std=c++0x")
	endif()
endif()
if(GLAM_CXX11_FLAGS)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLAM_CXX11_FLAGS}")
endif()

check_cxx_source_runs("
#include <tuple>
struct Test { };
int main() {
	std::tuple<int, float, Test> t = std::make_tuple(1, 0.5f, Test());
	return 0;
}
" HAS_TUPLE)
if (!HAS_TUPLE)
	message(FATAL_ERROR "Your compiler doesn't support C++11 tuples.")
endif()

check_cxx_source_runs("
int main() {
	static_assert(true, \"true is true\");
	return 0;
}
" HAS_STATIC_ASSERT)
if (!HAS_STATIC_ASSERT)
	message(FATAL_ERROR "Your compiler doesn't support C++11 static assert.")
endif()

check_cxx_source_runs("
#include <cmath>
int main() {
	std::round(0.5);
	return 0;
}
" HAS_CMATH)
if (!HAS_CMATH)
	message(FATAL_ERROR "Your compiler doesn't support C++11 <cmath>.")
endif()

check_cxx_source_runs("
void test(int i) { }
template <typename... Args>
void test(int i, Args... args) {
	test(args...);
}
int main() {
	test(1, 2, 3, 4, 5);
	return 0;
}
" HAS_VARIADIC_TEMPLATES)
if (!HAS_VARIADIC_TEMPLATES)
	message(FATAL_ERROR "Your compiler doesn't support C++11 variadic templates.")
endif()

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/glam/config.h"
)

enable_testing()

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(glam)
add_subdirectory(test)
add_subdirectory(example)
