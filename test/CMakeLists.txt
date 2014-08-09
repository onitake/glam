include_directories("${PROJECT_SOURCE_DIR}")

if(CXXTEST_FOUND)
	CXXTEST_ADD_TEST(math cxxtest_math.cpp math.cpp)
	CXXTEST_ADD_TEST(vector cxxtest_vector.cpp vector.cpp)
	CXXTEST_ADD_TEST(matrix cxxtest_matrix.cpp matrix.cpp)
endif()
