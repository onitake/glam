/* Copyright (c) 2014, Gregor Riepl <onitake@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <glam/vector.h>
#include <glam/matrix.h>
#include <glam/config.h>

#ifdef HAS_CXXTEST
#include <cxxtest/TestSuite.h>
#else
#include <iostream>
struct AssertException {
	unsigned long line;
	AssertException(unsigned long l) : line(l) { }
};
#define TS_ASSERT(a) do { if (!a) throw AssertException(__LINE__); };
#define TS_ASSERT_EQUALS(a, b) do { if (a != b) throw AssertException(__LINE__); };
#define TS_ASSERT_DELTA(a, b, d) do { if (a - b > d && b - a > d) throw AssertException(__LINE__); };
#endif

#ifdef HAS_CXXTEST
class VectorTest : public CxxTest::TestSuite {
#else
class VectorTest {
	VectorTest() {
		testVec4ConstructorAndAccessor();
		testVec5Equal();
		testVec4Multiply();
		testVec4Add();
		testIVec4ConstructorAndAccessor();
		testIVec4Cast();
		testVec4Dot();
		testDVec10IteratorConstructor();
		testDVec10Length();
		testDVec10Normalize();
	}
#endif
public:
	void testVec4ConstructorAndAccessor() {
		glam::vec4 v(1, 12.345, 0.5, 30);
		TS_ASSERT_EQUALS(v[0], 1.0f);
		TS_ASSERT_EQUALS(v[1], 12.345f);
		TS_ASSERT_EQUALS(v[2], 0.5f);
		TS_ASSERT_EQUALS(v[3], 30.0f);
	}
	void testVec4Assign() {
		glam::vec4 v(1, 12.345, 0.5, 30);
		glam::vec4 u = v;
		TS_ASSERT_EQUALS(u[0], 1.0f);
		TS_ASSERT_EQUALS(u[1], 12.345f);
		TS_ASSERT_EQUALS(u[2], 0.5f);
		TS_ASSERT_EQUALS(u[3], 30.0f);
	}
	void testVec5Equal() {
		glam::Vector<float, 5> v(1, 2, 3.5, 100, 0.00001);
		TS_ASSERT(v == v);
		glam::Vector<float, 5> v2(1, 2, 3.5, 100, 0.00001);
		TS_ASSERT(v == v2);
		glam::Vector<float, 5> w(0, 2, 3.5, 200, 0.0001);
		TS_ASSERT(v != w);
		glam::Vector<float, 5> w1(1, 2, 3.5, 100, 0.00001);
		TS_ASSERT(glam::all(glam::equal(v, w1)));
		glam::Vector<float, 5> w2(1, 2, 3.5, 100, 0.0001);
		TS_ASSERT(!glam::all(glam::equal(v, w2)));
		TS_ASSERT(glam::any(glam::equal(v, w2)));
		glam::Vector<float, 5> w3(0, 0, 0, 0, 0);
		TS_ASSERT(!glam::all(glam::equal(v, w3)));
		TS_ASSERT(!glam::any(glam::equal(v, w3)));
		TS_ASSERT(glam::all(glam::notEqual(v, w3)));
		TS_ASSERT(glam::any(glam::notEqual(v, w3)));
		glam::Vector<float, 5> w4(0.0001, 100, 3.5, 2, 1);
		TS_ASSERT(!glam::all(glam::equal(v, w4)));
		TS_ASSERT(glam::any(glam::equal(v, w4)));
		TS_ASSERT(glam::any(glam::notEqual(v, w4)));
	}
	void testVec4Multiply() {
		glam::vec4 v(1, 2, 3.5, 100);
		glam::vec4 w(0.5, 80, 20, 0.1);
		glam::vec4 u = v * w;
		TS_ASSERT_EQUALS(u[0], 0.5f);
		TS_ASSERT_EQUALS(u[1], 160.0f);
		TS_ASSERT_EQUALS(u[2], 70.0f);
		TS_ASSERT_EQUALS(u[3], 10.0f);
	}
	void testVec4Add() {
		glam::vec4 v(1, 2, 3.5, 100);
		glam::vec4 w(0.5, 80, 20, 0.1);
		glam::vec4 u = v + w;
		TS_ASSERT_EQUALS(u[0], 1.5f);
		TS_ASSERT_EQUALS(u[1], 82.0f);
		TS_ASSERT_EQUALS(u[2], 23.5f);
		TS_ASSERT_EQUALS(u[3], 100.1f);
	}
	void testIVec4ConstructorAndAccessor() {
		glam::ivec4 v(1, 12.345, 0.5, 30);
		TS_ASSERT_EQUALS(v[0], 1);
		TS_ASSERT_EQUALS(v[1], 12);
		TS_ASSERT_EQUALS(v[2], 0);
		TS_ASSERT_EQUALS(v[3], 30);
	}
	void testIVec4Cast() {
		glam::ivec4 v(1, 12.345, 0.5, 30);
		glam::vec4 w(v);
		TS_ASSERT_EQUALS(w[0], 1);
		TS_ASSERT_EQUALS(w[1], 12);
		TS_ASSERT_EQUALS(w[2], 0);
		TS_ASSERT_EQUALS(w[3], 30);
	}
	void testVec4Dot() {
		glam::vec4 v(1, 2, 3.5, 100);
		glam::vec4 w(0.5, 80, 20, 0.1);
		float u = glam::dot(v, w);
		TS_ASSERT_EQUALS(u, 240.5f);
	}
	void testDVec10IteratorConstructor() {
		double values[] = { 10, 100.1, 42, 9.856, 19, 0, 37.3, 90, 101.85, 0.00007 };
		glam::Vector<double, 10> v(std::make_pair(values, &values[10]));
		TS_ASSERT_EQUALS(v[0], 10.0);
		TS_ASSERT_EQUALS(v[1], 100.1);
		TS_ASSERT_EQUALS(v[2], 42.0);
		TS_ASSERT_EQUALS(v[3], 9.856);
		TS_ASSERT_EQUALS(v[4], 19.0);
		TS_ASSERT_EQUALS(v[5], 0.0);
		TS_ASSERT_EQUALS(v[6], 37.3);
		TS_ASSERT_EQUALS(v[7], 90.0);
		TS_ASSERT_EQUALS(v[8], 101.85);
		TS_ASSERT_EQUALS(v[9], 0.00007);
	}
	void testDVec10Length() {
		double values[] = { 10, 100.1, 42, 9.856, 19, 0, 37.3, 90, 101.85, 0.00007 };
		glam::Vector<double, 10> v(std::make_pair(values, &values[10]));
		double l = glam::length(v);
		TS_ASSERT_DELTA(l, 179.4627070898154, 1e-13);
	}
	void testDVec10Normalize() {
		double values[] = { 10, 100.1, 42, 9.856, 19, 0, 37.3, 90, 101.85, 0.00007 };
		glam::Vector<double, 10> v(std::make_pair(values, &values[10]));
		glam::Vector<double, 10> u = glam::normalize(v);
		TS_ASSERT_DELTA(u[0], 0.055721883182088, 1e-13);
		TS_ASSERT_DELTA(u[1], 0.5577760506527, 1e-13);
		TS_ASSERT_DELTA(u[2], 0.23403190936477, 1e-13);
		TS_ASSERT_DELTA(u[3], 0.054919488064266, 1e-13);
		TS_ASSERT_DELTA(u[4], 0.10587157804596, 1e-13);
		TS_ASSERT_DELTA(u[5], 0.0, 1e-13);
		TS_ASSERT_DELTA(u[6], 0.20784262426918, 1e-13);
		TS_ASSERT_DELTA(u[7], 0.50149694863879, 1e-13);
		TS_ASSERT_DELTA(u[8], 0.56752738020956, 1e-13);
		TS_ASSERT_DELTA(u[9], 0.90053182274617e-7, 1e-6);
	}
};

#ifndef HAS_CXXTEST
int main(int argc, char **argv) {
	try {
		VectorTest();
	} catch (AssertException &e) {
		std::cerr << "Assertion failed on line " << e.line << std::endl;
		return 1;
	} catch () {
		return 2;
	}
	return 0;
}
#endif
