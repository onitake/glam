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

#include <cxxtest/TestSuite.h>

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
const float M5[3*3] = {
	1, 2, 3,
	4, 5, 6,
	7, 8, 10,
};
const float M5D = -3;
const float M5I[3*3] = {
	-0.66666666666666f, -1.333333333333333f, 1.0f,
	-0.66666666666666f, 3.666666666666666f, -2.0f,
	1.0f, -2.0f, 1.0f,
};
const int M6[4*4] = {
	100, 33, 71, 19,
	21, 69, 98, 3,
	45, 28, 13, 74,
	62, 7, 19, 48,
};
const int V6[4] = {
	10, 20, 30, 40,
};
const int P6[4] = {
	5250, 2830, 3820, 4390,
};
const float M7A[4*2] = {
	0.3f, 0.5f,
	0.7f, 1.1f,
	1.3f, 1.7f,
	1.9f, 2.3f,
};
const float M7B[3*4] = {
	2.7f, 3.1f, 3.7f, 4.1f,
	4.3f, 4.7f, 5.1f, 5.3f,
	5.7f, 5.9f, 6.1f, 6.7f,
};
const float P7[3*2] = {
	15.58f, 20.48f,
	21.28f, 28.18f,
	26.5f, 35.12f,
};

class MatrixTest : public CxxTest::TestSuite {
public:
	void testMat2ConstructorAndAccessor() {
		glam::mat2 m2a(std::make_pair(M2, &M2[4]));
		TS_ASSERT_EQUALS(m2a[0][0], M2[0]);
		TS_ASSERT_EQUALS(m2a[0][1], M2[1]);
		TS_ASSERT_EQUALS(m2a[1][0], M2[2]);
		TS_ASSERT_EQUALS(m2a[1][1], M2[3]);
		glam::mat2 m2b(M2);
		TS_ASSERT_EQUALS(m2b[0][0], M2[0]);
		TS_ASSERT_EQUALS(m2b[0][1], M2[1]);
		TS_ASSERT_EQUALS(m2b[1][0], M2[2]);
		TS_ASSERT_EQUALS(m2b[1][1], M2[3]);
	}
	void testMat4IdentityAndZero() {
		glam::mat4 iden(1.0f);
		for (unsigned int j = 0; j < 4; j++) {
			for (unsigned int i = 0; i < 4; i++) {
				if (i == j) {
					TS_ASSERT_EQUALS(iden[j][i], 1.0f);
				} else {
					TS_ASSERT_EQUALS(iden[j][i], 0.0f);
				}
			}
		}
		glam::mat4 zero;
		for (unsigned int j = 0; j < 4; j++) {
			for (unsigned int i = 0; i < 4; i++) {
				TS_ASSERT_EQUALS(zero[j][i], 0.0f);
			}
		}
	}
	void testMat2Assign() {
		glam::mat2 m2(std::make_pair(M2, &M2[4]));
		glam::mat2 m2a = m2;
		TS_ASSERT_EQUALS(m2a[0][0], M2[0]);
		TS_ASSERT_EQUALS(m2a[0][1], M2[1]);
		TS_ASSERT_EQUALS(m2a[1][0], M2[2]);
		TS_ASSERT_EQUALS(m2a[1][1], M2[3]);
	}
	void testMat2Transpose() {
		glam::mat2 m2(std::make_pair(M2, &M2[4]));
		glam::mat2 m2t = glam::transpose(m2);
		TS_ASSERT_EQUALS(m2t[0][0], M2[0]);
		TS_ASSERT_EQUALS(m2t[0][1], M2[2]);
		TS_ASSERT_EQUALS(m2t[1][0], M2[1]);
		TS_ASSERT_EQUALS(m2t[1][1], M2[3]);
	}
	void testMat3DeterminantAndInverse() {
		glam::mat3 m5(std::make_pair(M5, &M5[9]));
		float d = glam::determinant(m5);
		TS_ASSERT_DELTA(d, M5D, 1e-6f);
		glam::mat3 m5i = glam::inverse(m5);
		TS_ASSERT_DELTA(m5i[0][0], M5I[0], 1e-6f);
		TS_ASSERT_DELTA(m5i[0][1], M5I[1], 1e-6f);
		TS_ASSERT_DELTA(m5i[0][2], M5I[2], 1e-6f);
		TS_ASSERT_DELTA(m5i[1][0], M5I[3], 1e-6f);
		TS_ASSERT_DELTA(m5i[1][1], M5I[4], 1e-6f);
		TS_ASSERT_DELTA(m5i[1][2], M5I[5], 1e-6f);
		TS_ASSERT_DELTA(m5i[2][0], M5I[6], 1e-6f);
		TS_ASSERT_DELTA(m5i[2][1], M5I[7], 1e-6f);
		TS_ASSERT_DELTA(m5i[2][2], M5I[8], 1e-6f);
	}
	void testIMat4Vector4Product() {
		glam::Matrix<int, 4, 4> m6(std::make_pair(M6, &M6[16]));
		glam::Vector<int, 4> v6(std::make_pair(V6, &V6[4]));
		glam::Vector<int, 4> p6 = m6 * v6;
		TS_ASSERT_EQUALS(p6[0], P6[0]);
		TS_ASSERT_EQUALS(p6[1], P6[1]);
		TS_ASSERT_EQUALS(p6[2], P6[2]);
		TS_ASSERT_EQUALS(p6[3], P6[3]);
	}
	void testMat42Mat23Product() {
		glam::mat4x2 m7a(std::make_pair(M7A, &M7A[8]));
		glam::mat3x4 m7b(std::make_pair(M7B, &M7B[12]));
		glam::mat3x2 p7 = m7a * m7b;
		TS_ASSERT_DELTA(p7[0][0], P7[0], 1e-6f);
		TS_ASSERT_DELTA(p7[0][1], P7[1], 1e-6f);
		TS_ASSERT_DELTA(p7[1][0], P7[2], 1e-6f);
		TS_ASSERT_DELTA(p7[1][1], P7[3], 1e-6f);
		TS_ASSERT_DELTA(p7[2][0], P7[4], 1e-6f);
		TS_ASSERT_DELTA(p7[2][1], P7[5], 1e-6f);
	}
};
