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
		long double d3 = 180.0L;
		long double r3 = glam::radians(d3);
		TS_ASSERT_DELTA(r3, 3.141592653589793238462643383279502884L, 1e-24);
		std::complex<double> d4 = 180.0;
		std::complex<double> r4 = glam::radians(d4);
		TS_ASSERT_DELTA(r4.real(), 3.141592653589793, 1e-12);
		TS_ASSERT_DELTA(r4.imag(), 0.0, 1e-12);
	}
};
