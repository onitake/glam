include_directories("${PROJECT_SOURCE_DIR}")

add_executable(vecmath vecmath.cpp)
add_executable(glslvec glslvec.cpp)
add_executable(linalgperf linalgperf.cpp)
add_executable(dft dft.cpp)

add_test(VectorMathTest vecmath)
add_test(GLSLMathTest glslvec)
add_test(LinearAlgebraPerformance linalgperf)
add_test(DFT dft)

