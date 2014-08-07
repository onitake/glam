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

#include <complex>
#include <glam/math.h>
#include <glam/config.h>

#include <cxxtest/TestSuite.h>

class MathTest : public CxxTest::TestSuite {
public:
	void testConstants() {
		TS_ASSERT_EQUALS(glam::CONSTANTS<int>::PI, 3);
		TS_ASSERT_DELTA(glam::CONSTANTS<float>::PI, 3.141592653589793f, 1e-6);
		TS_ASSERT_DELTA(glam::CONSTANTS<double>::PI, 3.141592653589793, 1e-12);
#ifdef GLAM_HAS_LONG_DOUBLE
		TS_ASSERT_DELTA(glam::CONSTANTS<long double>::PI, 3.141592653589793238462643383279502884L, 1e-24);
#endif
#ifndef GLAM_DEBUG
		TS_ASSERT_DELTA(glam::CONSTANTS<std::complex<double> >::PI.real(), 3.141592653589793, 1e-12);
		TS_ASSERT_DELTA(glam::CONSTANTS<std::complex<double> >::PI.imag(), 0.0, 1e-12);
#endif
		TS_ASSERT_EQUALS(glam::CONSTANTS<int>::E, 2);
		TS_ASSERT_DELTA(glam::CONSTANTS<float>::E, 2.7182818284590452354f, 1e-6);
		TS_ASSERT_DELTA(glam::CONSTANTS<double>::E, 2.7182818284590452354, 1e-12);
#ifdef GLAM_HAS_LONG_DOUBLE
		TS_ASSERT_DELTA(glam::CONSTANTS<long double>::E, 2.718281828459045235360287471352662498L, 1e-24);
#endif
#ifndef GLAM_DEBUG
		TS_ASSERT_DELTA(glam::CONSTANTS<std::complex<double> >::E.real(), 2.7182818284590452354, 1e-12);
		TS_ASSERT_DELTA(glam::CONSTANTS<std::complex<double> >::E.imag(), 0.0, 1e-12);
#endif
	}
	void testRadians() {
		int d0 = 180.0f;
		int r0 = glam::radians(d0);
		TS_ASSERT_EQUALS(r0, 3);
		float d1 = 180.0f;
		float r1 = glam::radians(d1);
		TS_ASSERT_DELTA(r1, 3.141592653589793f, 1e-6);
		double d2 = 180.0;
		double r2 = glam::radians(d2);
		TS_ASSERT_DELTA(r2, 3.141592653589793, 1e-12);
#ifdef GLAM_HAS_LONG_DOUBLE
		long double d3 = 180.0L;
		long double r3 = glam::radians(d3);
		TS_ASSERT_DELTA(r3, 3.141592653589793238462643383279502884L, 1e-24);
#endif
#ifndef GLAM_DEBUG
		std::complex<double> d4 = 180.0;
		std::complex<double> r4 = glam::radians(d4);
		TS_ASSERT_DELTA(r4.real(), 3.141592653589793, 1e-12);
		TS_ASSERT_DELTA(r4.imag(), 0.0, 1e-12);
#endif
	}
	void testBits() {
		int i = glam::floatBitsToInt(1.0f);
		TS_ASSERT_EQUALS(i, 0x3f800000);
		unsigned int u = glam::floatBitsToUInt(1.0f);
		TS_ASSERT_EQUALS(u, 0x3f800000);
		float f1 = glam::intBitsToFloat(0x3f800000);
		TS_ASSERT_EQUALS(f1, 1.0f);
		float f2 = glam::uintBitsToFloat(0x3f800000);
		TS_ASSERT_EQUALS(f2, 1.0f);
	}
};
