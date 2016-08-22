/*
 * GLAM - GLSL Linear Algebra Math Library
 * 
 * Copyright (c) 2014, Gregor Riepl <onitake@gmail.com>
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

#include <cxxtest/TestSuite.h>
#include <string>
#include <type_traits>
#include <glam/vector.h>
#include <glam/matrix.h>
#include <glam/config.h>
#include <cstdlib>
#include <cmath>

class VectorTest : public CxxTest::TestSuite {
public:
	void testSize() {
		TS_TRACE(std::string("sizeof(float) is ") + std::to_string(sizeof(float)));
		TS_TRACE(std::string("sizeof(vec2) is ") + std::to_string(sizeof(glam::vec2)));
		TS_TRACE(std::string("sizeof(vec4) is ") + std::to_string(sizeof(glam::vec4)));
		TS_TRACE(std::string("sizeof(bool) is ") + std::to_string(sizeof(bool)));
		TS_TRACE(std::string("sizeof(bvec4) is ") + std::to_string(sizeof(glam::bvec4)));
		TS_TRACE(std::string("sizeof(int) is ") + std::to_string(sizeof(int)));
		TS_TRACE(std::string("sizeof(ivec2) is ") + std::to_string(sizeof(glam::ivec2)));
		TS_TRACE(std::string("sizeof(double) is ") + std::to_string(sizeof(double)));
		TS_TRACE(std::string("sizeof(dvec3) is ") + std::to_string(sizeof(glam::dvec3)));
		// A vectorized implementation may pad to up to 4 vector elements to improve SIMD performance, but not more
		TS_ASSERT_LESS_THAN_EQUALS(sizeof(glam::Vector<float, 4>), sizeof(float) * 4);
		TS_ASSERT_LESS_THAN_EQUALS(sizeof(glam::Vector<double, 2>), sizeof(double) * 4);
		TS_ASSERT_LESS_THAN_EQUALS(sizeof(glam::Vector<bool, 4>), sizeof(bool) * 4);
		TS_ASSERT_LESS_THAN_EQUALS(sizeof(glam::Vector<int, 256>), sizeof(int) * 256);
	}
	void testFeatures() {
		// is_pod requires non-trivial copy/move constructors and non-trivial copy/move assignment operators
		// since we already provide a const internal pointer accessor, we don't enforce this requirementi
		TS_ASSERT(std::is_standard_layout<glam::vec2>::value);
		TS_ASSERT(std::is_standard_layout<glam::vec3>::value);
		TS_ASSERT(std::is_standard_layout<glam::vec4>::value);
		TS_ASSERT(std::is_standard_layout<glam::ivec2>::value);
		TS_ASSERT(std::is_standard_layout<glam::ivec3>::value);
		TS_ASSERT(std::is_standard_layout<glam::ivec4>::value);
		TS_ASSERT(std::is_standard_layout<glam::dvec2>::value);
		TS_ASSERT(std::is_standard_layout<glam::dvec3>::value);
		TS_ASSERT(std::is_standard_layout<glam::dvec4>::value);
		TS_ASSERT(std::is_standard_layout<glam::bvec2>::value);
		TS_ASSERT(std::is_standard_layout<glam::bvec3>::value);
		TS_ASSERT(std::is_standard_layout<glam::bvec4>::value);
	}
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
		TS_ASSERT(!(v != v));
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
	void testBVec4Logic() {
		glam::bvec4 v(true, false, true, false);
		glam::bvec4 u = not(v);
		TS_ASSERT(!u[0]);
		TS_ASSERT(u[1]);
		TS_ASSERT(!u[2]);
		TS_ASSERT(u[3]);
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
		glam::vec4 u2(v);
		glam::vec4 u1 = v + w;
		TS_ASSERT_EQUALS(u1[0], 1.5f);
		TS_ASSERT_EQUALS(u1[1], 82.0f);
		TS_ASSERT_EQUALS(u1[2], 23.5f);
		TS_ASSERT_EQUALS(u1[3], 100.1f);
		u2 += w;
		TS_ASSERT_EQUALS(u2[0], 1.5f);
		TS_ASSERT_EQUALS(u2[1], 82.0f);
		TS_ASSERT_EQUALS(u2[2], 23.5f);
		TS_ASSERT_EQUALS(u2[3], 100.1f);
	}
	void testIVec10Subtract() {
		glam::Vector<int, 10> v(1, 2, 3, 100, 4, 5, 6, 200, 7, 8);
		glam::Vector<int, 10> w(1000, 9, 8, 7, 800, 6, 5, 4, 600, 3);
		glam::Vector<int, 10> u2(v);
		glam::Vector<int, 10> u1 = v - w;
		TS_ASSERT_EQUALS(u1[0], -999);
		TS_ASSERT_EQUALS(u1[1], -7);
		TS_ASSERT_EQUALS(u1[2], -5);
		TS_ASSERT_EQUALS(u1[3], 93);
		TS_ASSERT_EQUALS(u1[4], -796);
		TS_ASSERT_EQUALS(u1[5], -1);
		TS_ASSERT_EQUALS(u1[6], 1);
		TS_ASSERT_EQUALS(u1[7], 196);
		TS_ASSERT_EQUALS(u1[8], -593);
		TS_ASSERT_EQUALS(u1[9], 5);
		u2 -= w;
		TS_ASSERT_EQUALS(u2[0], -999);
		TS_ASSERT_EQUALS(u2[1], -7);
		TS_ASSERT_EQUALS(u2[2], -5);
		TS_ASSERT_EQUALS(u2[3], 93);
		TS_ASSERT_EQUALS(u2[4], -796);
		TS_ASSERT_EQUALS(u2[5], -1);
		TS_ASSERT_EQUALS(u2[6], 1);
		TS_ASSERT_EQUALS(u2[7], 196);
		TS_ASSERT_EQUALS(u2[8], -593);
		TS_ASSERT_EQUALS(u2[9], 5);
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
	void testVecDot() {
		glam::vec4 v1(1, 2, 3.5, 100);
		glam::vec4 w1(0.5, 80, 20, 0.1);
		float u1 = glam::dot(v1, w1);
		TS_ASSERT_EQUALS(u1, 240.5f);
		glam::Vector<double, 7> v2(1, 2, 3, 4, 5, 6, 0);
		glam::Vector<double, 7> w2(0, 6, 5, 4, 3, 2, 1);
		TS_ASSERT_EQUALS(glam::dot(v2, w2), 70);
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
	void testVec2Permute() {
		glam::vec2 v(1, 2);
		TS_ASSERT_EQUALS(v.x, 1.0f);
		TS_ASSERT_EQUALS(v.y, 2.0f);
		glam::vec2 u1 = v.xy;
		TS_ASSERT_EQUALS(u1[0], 1.0f);
		TS_ASSERT_EQUALS(u1[1], 2.0f);
		glam::vec2 u2 = v.yx;
		TS_ASSERT_EQUALS(u2[0], 2.0f);
		TS_ASSERT_EQUALS(u2[1], 1.0f);
		glam::vec2 u3 = v.xx;
		TS_ASSERT_EQUALS(u3[0], 1.0f);
		TS_ASSERT_EQUALS(u3[1], 1.0f);
		glam::vec2 u4 = v.yy;
		TS_ASSERT_EQUALS(u4[0], 2.0f);
		TS_ASSERT_EQUALS(u4[1], 2.0f);
		glam::vec2 u5 = v;
		u5.y = 5.0f;
		TS_ASSERT_EQUALS(u5[0], 1.0f);
		TS_ASSERT_EQUALS(u5[1], 5.0f);
		glam::vec2 u6 = v;
		u6.yx = glam::vec2(3, 4);
		TS_ASSERT_EQUALS(u6[0], 4.0f);
		TS_ASSERT_EQUALS(u6[1], 3.0f);
		glam::vec2 u7 = v;
		glam::vec2 u8(10, 20);
		u7.yx = u8.yx;
		TS_ASSERT_EQUALS(u7[0], 10.0f);
		TS_ASSERT_EQUALS(u7[1], 20.0f);
	}
	void testVec3PermuteAlt() {
		glam::vec3 v;
		v.rgb = glam::vec3(0.5f, 0.0f, 1.0f);
		glam::vec3 u2 = v.bgr;
		TS_ASSERT_EQUALS(u2[0], 1.0f);
		TS_ASSERT_EQUALS(u2[1], 0.0f);
		TS_ASSERT_EQUALS(u2[2], 0.5f);
		glam::vec3 u3;
		u3.zyx = v.rgb;
		TS_ASSERT_EQUALS(u3[0], 1.0f);
		TS_ASSERT_EQUALS(u3[1], 0.0f);
		TS_ASSERT_EQUALS(u3[2], 0.5f);
		glam::vec3 u4;
		u4.stp = v.xyz;
		TS_ASSERT_EQUALS(u4[0], 0.5f);
		TS_ASSERT_EQUALS(u4[1], 0.0f);
		TS_ASSERT_EQUALS(u4[2], 1.0f);
		glam::vec3 u5;
		u5.p = v.x;
		TS_ASSERT_EQUALS(u5[0], 0.0f);
		TS_ASSERT_EQUALS(u5[1], 0.0f);
		TS_ASSERT_EQUALS(u5[2], 0.5f);
		glam::vec3 u6 = v;
		u6.tp = v.xx;
		TS_ASSERT_EQUALS(u6[0], 0.5f);
		TS_ASSERT_EQUALS(u6[1], 0.5f);
		TS_ASSERT_EQUALS(u6[2], 0.5f);
	}
	void testVec2PermuteMath() {
		glam::vec2 v(1, 2);
		glam::vec2 u1 = v;
		glam::vec2 u2(10, 20);
		u1.yx += u2;
		TS_ASSERT_EQUALS(u1[0], 21.0f);
		TS_ASSERT_EQUALS(u1[1], 12.0f);
	}
	void testVec3Bits() {
		glam::vec3 vf1(1.0f, -1.0f, 0.5f);
		glam::ivec3 vi1 = glam::floatBitsToInt(vf1);
		TS_ASSERT_EQUALS(vi1[0], 0x3f800000);
		TS_ASSERT_EQUALS(vi1[1], 0xbf800000);
		TS_ASSERT_EQUALS(vi1[2], 0x3f000000);
		glam::ivec3 vi2(0x3f800000, 0xbf800000, 0x3f000000);
		glam::vec3 vf2 = glam::intBitsToFloat(vi2);
		TS_ASSERT_EQUALS(vf2[0], 1.0f);
		TS_ASSERT_EQUALS(vf2[1], -1.0f);
		TS_ASSERT_EQUALS(vf2[2], 0.5f);
	}
	void testPackNorm() {
		glam::vec2 vsf1(-1.0f, 1.0f);
		unsigned int up1 = glam::packSnorm2x16(vsf1);
		// Should be 8000, precision? Rounding error?
		TS_ASSERT_EQUALS(up1, 0x7fff8001);
		unsigned int up2 = 0x7fff8000;
		glam::vec2 vsf2 = glam::unpackSnorm2x16(up2);
		TS_ASSERT_DELTA(vsf2[0], -1.0f, 1e-6);
		TS_ASSERT_DELTA(vsf2[1], 1.0f, 1e-6);
		glam::vec2 vuf1(0.0f, 1.0f);
		unsigned int up3 = glam::packUnorm2x16(vuf1);
		TS_ASSERT_EQUALS(up3, 0xffff0000);
		unsigned int up4 = 0xffff0000;
		glam::vec2 vuf2 = glam::unpackUnorm2x16(up4);
		TS_ASSERT_DELTA(vuf2[0], 0.0f, 1e-6);
		TS_ASSERT_DELTA(vuf2[1], 1.0f, 1e-6);
	}
	void testPackHalf() {
		glam::vec2 vf1(10.0f, -0.125f);
		unsigned int p1 = glam::packHalf2x16(vf1);
		TS_ASSERT_EQUALS(p1, 0xb0004900);
		unsigned int p2 = 0xb0004900;
		glam::vec2 vf2 = glam::unpackHalf2x16(p2);
		TS_ASSERT_DELTA(vf2[0], 10.0f, 1e-4);
		TS_ASSERT_DELTA(vf2[1], -0.125f, 1e-4);
		glam::vec2 vf3(768.0f, 5.779e-41f);
		unsigned int p3 = glam::packHalf2x16(vf3);
		TS_ASSERT_EQUALS(p3, 0x00056200);
		unsigned int p4 = 0x00056200;
		glam::vec2 vf4 = glam::unpackHalf2x16(p4);
		TS_ASSERT_DELTA(vf4[0], 768.0f, 1e-4);
		TS_ASSERT_DELTA(vf4[1], 5.779e-41f, 1e-4);
		glam::vec2 vf5(-INFINITY, NAN);
		unsigned int p5 = glam::packHalf2x16(vf5);
		//std::printf("vf5=(%g,%g)=(0x%08x,0x%08x) p5=0x%08x\n", vf5.x, vf5.y, *(uint32_t *) &vf5.x, *(uint32_t *) &vf5.y, p5);
		TS_ASSERT((p5 & 0xffff) == 0xfc00 && (p5 & 0x7c000000) == 0x7c000000 && (p5 & 0x03ff0000) != 0x00000000);
		unsigned int p6 = 0x7e00fc00;
		glam::vec2 vf6 = glam::unpackHalf2x16(p6);
		//std::printf("vf6=(%g,%g)=(0x%08x,0x%08x) p6=0x%08x\n", vf6.x, vf6.y, *(uint32_t *) &vf6.x, *(uint32_t *) &vf6.y, p6);
#ifndef GLAM_HAS_BROKEN_ISINF
		// this will fail with -ffast-math, according to http://stackoverflow.com/questions/22931147/22931368#22931368
		TS_ASSERT(std::isinf(vf6[0]) && std::signbit(vf6[0]) && std::isnan(vf6[1]));
#else
		TS_ASSERT(*(uint32_t *) &vf6[0] == 0xff800000 && *(uint32_t *) &vf6[1] == 0x7fc00000);
#endif
		for (int i = 0; i <= 15; i++) {
			glam::vec2 vfi(glam::pow(1.5f, float(i)), glam::pow(2.0f, -float(i)));
			unsigned int pi = glam::packHalf2x16(vfi);
			glam::vec2 vfj = glam::unpackHalf2x16(pi);
			//std::printf("vfi=(%g,%g)=(0x%08x,0x%08x) pi=0x%08x vfj=(%g,%g)=(0x%08x,0x%08x)\n", vfi.x, vfi.y, *(uint32_t *) &vfi.x, *(uint32_t *) &vfi.y, pi, vfj.x, vfj.y, *(uint32_t *) &vfj.x, *(uint32_t *) &vfj.y);
			TS_ASSERT_DELTA(vfi[0], vfj[0], 1e0);
			TS_ASSERT_DELTA(vfi[1], vfj[1], 1e-4);
		}
	}
};
