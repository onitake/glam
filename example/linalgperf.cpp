/* Copyright (c) 2012-2014, Gregor Riepl <onitake@gmail.com>
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

#include <iostream>
#include <cassert>
#include <cmath>
#include <stdint.h>
#include <glam/matrix.h>

using glam::determinant;
using glam::inverse;
using glam::transpose;
typedef glam::Matrix<float, 2, 2> Matrix2f;
typedef glam::Matrix<float, 5, 5> Matrix5f;
typedef glam::Matrix<float, 6, 6> Matrix6f;
typedef glam::Matrix<float, 10, 10> Matrix10f;
typedef glam::Matrix<double, 10, 10> Matrix10d;

// Transposed notation
const float M1[6*6] = {
	1, 2, 3, 4, 5, 6,
	7, 8, 9, 10, 11, 12,
	13, 14, 15, 16, 17, 18,
	19, 20, 21, 22, 23, 24,
	25, 26, 27, 28, 29, 30,
	31, 32, 33, 34, 35, 36,
};
const float M2[2*2] = {
	1, 2,
	3, 4,
};
const float M3[5*5] = {
	152, 218, 190, 27, 46,
	222, 225, 181, 110, 120,
	126, 230, 28, 29, 19,
	248, 69, 186, 160, 224,
	154, 173, 197, 49, 43,
};
// The 0 allows testing if permutations work
const double M4[10*10] = {
	0, 94, 23, 152, 91, 248, 95, 250, 186, 165, 
	69, 61, 56, 182, 247, 14, 253, 146, 32, 92, 
	215, 247, 75, 192, 117, 117, 111, 118, 7, 38, 
	1, 55, 49, 162, 13, 255, 135, 142, 145, 247, 
	93, 15, 121, 46, 103, 138, 118, 144, 67, 214, 
	96, 232, 114, 87, 229, 109, 180, 155, 188, 12, 
	13, 146, 36, 155, 184, 69, 59, 161, 232, 224, 
	81, 230, 54, 207, 125, 169, 1, 236, 135, 63, 
	74, 46, 184, 164, 253, 254, 42, 225, 229, 183, 
	214, 41, 224, 111, 125, 214, 22, 2, 92, 24, 
};

template Matrix2f glam::inverse(const Matrix2f &);
template Matrix5f glam::inverse(const Matrix5f &);
template Matrix6f glam::inverse(const Matrix6f &);
template Matrix10f glam::inverse(const Matrix10f &);
template Matrix10d glam::inverse(const Matrix10d &);

// Fuzzy comparison
// Adapted from http://randomascii.wordpress.com/2012/01/11/tricks-with-the-floating-point-format/
bool compare(float a, float b) {
	union Float {
		Float(float num = 0.0f) : f(num) { }
		// Portable extraction of components.
		bool negative() const { return (i >> 31) != 0; }
		int32_t fraction() const { return i & ((1 << 23) - 1); }
		int32_t exponent() const { return (i >> 23) & 0xff; }
		int32_t i;
		float f;
	};
	Float ua(a);
	Float ub(b);
	//std::cout << "a=" << (ua.negative() ? "-" : "") << ua.fraction() << "e" << (ua.exponent() - 127) << std::endl;
	//std::cout << "b=" << (ub.negative() ? "-" : "") << ub.fraction() << "e" << (ub.exponent() - 127) << std::endl;
	// Compare the exponents
	int32_t edif = std::abs(ua.exponent() - ub.exponent());
	//std::cout << "edif=" << edif << std::endl;
	if (edif == 0) {
		// Exponents are equal, allow for a small error in value
		int32_t fdif = std::abs(ua.fraction() - ub.fraction());
		//std::cout << "fdif=" << fdif << std::endl;
		if (fdif > 0x000100) {
			return false;
		}
		return true;
	}
	// Scale the relative allowable error according to the exponent of the smaller of a and b
	// Note that we do not subtract 127 from the exponent as per the definition of IEEE 754 single-precision
	// 128 and up: lowest allowable error
	// 1: maximum allowable error
	int32_t max = 128 - std::min(ua.exponent(), ub.exponent());
	if (max < 0) max = 0;
	//std::cout << "max=" << max << std::endl;
	if (edif > max) {
		return false;
	}
	return true;
}
bool compare(double a, double b) {
	union Double {
		Double(float num = 0.0f) : d(num) { }
		// Portable extraction of components.
		bool negative() const { return (i >> 63) != 0; }
		int64_t fraction() const { return i & ((1L << 52) - 1); }
		int64_t exponent() const { return (i >> 52) & 0x7ffL; }
		int64_t i;
		double d;
	};
	Double ua(a);
	Double ub(b);
	// Compare the exponents
	int64_t edif = std::abs(ua.exponent() - ub.exponent());
	//std::cout << "edif=" << edif << std::endl;
	if (edif == 0) {
		// Exponents are equal, allow for a small error in value
		int64_t fdif = std::abs(ua.fraction() - ub.fraction());
		//std::cout << "fdif=" << fdif << std::endl;
		if (fdif > 0x0000000010000L) {
			return false;
		}
		return true;
	}
	// Scale the relative allowable error according to the exponent of the smaller of a and b
	// Note that we do not subtract 1023 from the exponent as per the definition of IEEE 754 double-precision
	// 1024 and up: lowest allowable error
	// 1: maximum allowable error
	int64_t max = 1024L - std::min(ua.exponent(), ub.exponent());
	//std::cout << "max=" << max << std::endl;
	if (edif > max) {
		return false;
	}
	return true;
}

// Matrix equality, allowing for error
template <class T, unsigned int M, unsigned int N>
bool compare(const glam::Matrix <T, M, N> &a, const glam::Matrix <T, M, N> &b) {
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			if (!compare(a[i][j], b[i][j])) {
				return false;
			}
		}
	}
	return true;
}

int main(int argc, char **argv) {
	Matrix6f m1(M1);
	Matrix2f m2(M2);
	Matrix5f m3(M3);
	Matrix10d m4(M4);
	Matrix2f i2(1);
	Matrix5f i5(1);
	Matrix6f i6(1);
	Matrix10d i10(1);
	std::cout << m1 << std::endl;
	std::cout << transpose(m1) << std::endl;
	assert(m1 * i6 == m1);
	assert(i6 * m1 == m1);
	std::cout << m2 << std::endl;
	std::cout << inverse(m2) << std::endl;
	assert(m2 * inverse(m2) == i2);
	assert(inverse(m2) * m2 == i2);
	assert(determinant(m1) == 0);
	bool e = false;
	try {
		inverse(m1);
	} catch (glam::NonInvertibleMatrixException &ex) {
		e = true;
	}
	assert(e);
	std::cout << m3 << std::endl;
	std::cout << inverse(m3) << std::endl;
	std::cout << (m3 * inverse(m3)) << std::endl;
	assert(compare(m3 * inverse(m3), i5));
	assert(compare(inverse(m3) * m3, i5));
	assert(!compare(inverse(m2), m2));
	assert(!compare(inverse(m3), m3));
	//assert(compare(m3.inv(), m3.invr()));
	std::cout << "det(m3)[1] = " << determinant(m3) << std::endl;
	//std::cout << "det(m3)[2] = " << m3.detr() << std::endl;
	//assert(compare(m3.det(), m3.detr()));
	//Matrix10d m4i = m4.inv();
	//std::cout << m4i << std::endl;
	//std::cout << (m4 * m4i) << std::endl;
	//assert(compare(m4 * m4i, i10));
	/*std::vector <Matrix10d> lup = m4.lup();
	std::cout << lup[0] << std::endl;
	std::cout << lup[1] << std::endl;
	std::cout << lup[2] << std::endl;
	std::cout << m4 << std::endl;
	std::cout << (lup[0] * lup[1]) << std::endl;
	std::cout << (lup[1] * lup[0]) << std::endl;
	// A = P^-1LU, P^-1 = Pt
	assert(compare(lup[2].t() * lup[0] * lup[1], m4));
	assert(!compare(lup[1] * lup[0], m4));*/
	Matrix10d m4i = inverse(m4);
	std::cout << m4i << std::endl;
	assert(compare(m4 * m4i, i10));
	assert(!compare(m4, m4i));
	//Matrix10d m4i2 = m4.invr();
	//std::cout << m4i2 << std::endl;
	//assert(compare(m4 * m4i2, i10));
	return 0;
}

