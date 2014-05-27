include_directories("${PROJECT_SOURCE_DIR}")

if(CXXTEST_FOUND)
	CXXTEST_ADD_TEST(vector cxxtest_vector.cpp vector.cpp)
	CXXTEST_ADD_TEST(matrix cxxtest_matrix.cpp matrix.cpp)
else()
	add_executable(vector vector.cpp)
	add_test(VectorTest vector)
endif()
