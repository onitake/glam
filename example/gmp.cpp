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

// This example implements the test case from http://www.ams.org/notices/201410/rnoti-p1249.pdf,
// which produces incorrect results in Mathematica 8 and 9.
// Big integer arithmetics are provided by libgmp(xx).
// For self-verification purposes (and performance analysis), a run of 1000 calculations is performed.
// The result should be roughly 1.95124e+9762 and no errors should be reported.
// Note the matrix transposition before calculating the determinant, there seems to be an issue with
// column/row ordering. Additional verification needed.

#include <iostream>
#include <gmpxx.h>
#include <glam/vector.h>
#include <glam/matrix.h>

using namespace glam;
using std::cout;
using std::endl;
using std::make_pair;

typedef mpz_class BigInt;
typedef mpf_class BigFloat;
typedef Vector<unsigned long, 14> VectorULong14;
typedef Vector<mpz_class, 14> VectorBigInt14;
typedef Vector<mpf_class, 14> VectorBigFloat14;
typedef Matrix<mpz_class, 14, 14> MatrixBigInt14;
typedef Matrix<mpf_class, 14, 14> MatrixBigFloat14;

// gmp does not have a power function that takes two mpz arguments and no modulus, so we use the mpz^ulong variant instead
static inline VectorBigInt14 pow(const VectorBigInt14 &x, const VectorULong14 &y) {
	VectorBigInt14 ret;
	for (size_t c = 0; c < 14; c++) {
		mpz_pow_ui(ret[c].get_mpz_t(), x[c].get_mpz_t(), y[c]);
	}
	return ret;
}


static const unsigned long powersData[] = { 123, 152, 185, 220, 397, 449, 503, 563, 979, 1059, 1143, 1229, 1319, 1412 };
static const BigInt basicData[] = {
	-32, 6, 78, -93, 54, 65, 29, 93, 45, -13, 11, -72, -92, 47,
	69, 99, 38, 30, -22, -58, -58, -77, -99, 51, 72, 53, -22, -73,
	89, 11, 12, 89, 54, 3, 47, 16, 3, -99, -74, 49, -90, -30,
	-60, 57, 64, 22, -53, 57, 87, -64, -57, -2, 86, -46, 67, -60,
	-83, 47, -67, 13, -52, 19, 3, 48, 32, 48, 79, 17, -25, 99,
	-22, -42, -4, 48, 64, 77, -6, 84, 60, 71, -58, -22, -28, 9,
	-14, -48, -52, -73, 19, 76, -81, 97, 74, -81, -89, -48, -91, -86,
	-58, -65, -65, 93, 1, -57, 5, 75, 4, -32, 80, -40, -8, -70,
	85, 25, 19, 11, 81, -80, 98, 89, 69, 78, 70, -28, 32, 84,
	56, 50, 71, -97, -72, 22, 86, 63, 98, 27, 55, -85, -41, 55,
	-65, -70, 38, -49, -11, 93, -98, 34, -40, -28, -49, 88, 10, 19,
	-30, -3, -17, 61, 50, -85, 51, -98, -69, -22, 51, -30, 6, 69,
	-86, -90, 51, -25, 0, 67, -62, -94, -28, 22, -42, 74, 85, 11,
	-9, 31, -3, -4, -81, 58, -66, 19, -26, 94, 66, 32, 21, -84,
};
static const BigInt smallData[] = {
	528, 878, -357, -76, -358, 192, -234, 80, 967, 406, 33, -537, -373, 638,
	853, 910, 451, 26, 18, 233, -71, 748, 1, -944, -328, -311, -447, -33,
	-547, -306, -475, -778, -229, 620, -831, 684, -494, -533, -543, 829, -252, 932,
	-323, -260, 327, 505, -880, 955, 880, 332, 633, -827, 583, 28, -619, -335,
	393, 575, -84, 942, -955, -877, -135, 730, 891, 615, -443, 175, -418, -75,
	-916, -765, 237, -561, -346, 281, -249, -111, -907, 907, -635, 182, 764, -676,
	-11, -32, 647, -350, 550, 357, -427, -643, -586, -443, 904, -930, 994, -934,
	-976, 94, 505, 698, -958, -226, 737, 102, 129, -350, -745, 258, -543, 239,
	279, 254, -137, -532, 867, -820, 664, -242, 688, 700, -398, -808, -37, 210,
	-665, 276, 363, -507, -541, 513, 298, -82, 150, -878, -110, -399, -845, 665,
	906, -156, -808, -78, -962, -882, -552, -28, -501, 706, 751, -43, 30, 414,
	-277, 625, 332, -758, 646, 536, -1, 585, -298, 1, 660, -68, -704, -803,
	103, -8, 222, 346, 932, -237, -712, 207, 704, 800, 474, -553, 147, 564,
	-485, -566, -998, -545, 168, 877, -691, -986, -68, 120, 255, 421, -534, -805,
};

int main(int argc, char **argv) {
	bool set = false;
	BigInt alldet;
	for (int i = 0; i < 1000; i++) {
		MatrixBigInt14 basicMatrix = MatrixBigInt14(make_pair(basicData, basicData + 14 * 14));
		VectorULong14 powersVector = VectorULong14(make_pair(powersData, powersData + 14));
		MatrixBigInt14 powersMatrix = scalingMatrix(pow(VectorBigInt14(10), powersVector));
		MatrixBigInt14 smallMatrix = MatrixBigInt14(make_pair(smallData, smallData + 14 * 14));
		MatrixBigInt14 bigMatrix = transpose(basicMatrix * powersMatrix + smallMatrix);
		BigInt det = determinant(bigMatrix);
		if (set) {
			if (det != alldet) {
				cout << "Error, iteration " << i << " produced a different result: " << BigFloat(det) << endl;
			}
		} else {
			cout << "basicMatrix: " << basicMatrix << endl;
			cout << "powersVector: " << powersVector << endl;
			cout << "powersMatrix: " << MatrixBigFloat14(powersMatrix) << endl;
			cout << "smallMatrix: " << smallMatrix << endl;
			cout << "result: " << BigFloat(det) << endl;
			alldet = det;
			set = true;
		}
	}
	return 0;
}
