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

#ifndef GLAM_MATRIX_H
#define GLAM_MATRIX_H

#include <cstring>
#include <cmath>
#include <ostream>
#include <utility>
#include <tuple>
#include <glam/vector.h>
#include <glam/exception.h>
#include <glam/config.h>

namespace glam {

// Forward declarations
template <class T, unsigned int M, unsigned int N>
class Matrix;

// M=height N=width i=row j=column
template <class T, unsigned int M, unsigned int N = M>
class Matrix {
private:
	// The matrix is sorted in column-major order, as in OpenGL:
	// | a00 a01 a02 | => [ a00 a10 ] [ a01 a11 ] [ a02 a12 ] (index = M * j + i)
	// | a10 a11 a12 |
	Vector<T, M> _a[N];

	// Single element access by array index
	T &element(unsigned int index);
	const T &element(unsigned int index) const;
	// LU(P) decomposition
	// LOWER is a lower triangular matrix
	// UPPER is an upper triangular matrix
	// PERMUTATION is a permutation matrix
	// SWAPS is the number of swaps performed to obtain PERMUTATION, with respect to IDENTITY
	// The following equation must hold: PA = LU (where A is *this)
	enum TupleElement { LOWER = 0, UPPER = 1, PERMUTATION = 2, SWAPS = 3 };
	std::tuple<Matrix<T, M, N>, Matrix<T, M, N>, Matrix<T, M, N>, unsigned int> lu() const;

	// Supporting initializers to avoid trouble from delegate constructors
	template <unsigned int X>
	void initialize();
	template <unsigned int X>
	void initialize(const T &v);
	template <unsigned int X, typename ForwardIterator>
	void initialize(std::pair<ForwardIterator, ForwardIterator> iterator);
	template <unsigned int X>
	void initialize(const T *v);
	template <unsigned int X, class U, unsigned int P, unsigned int R>
	void initialize(const Matrix<U, P, R> &m);
	template <unsigned int X, typename... Args, unsigned int P>
	void initialize(const Vector<T, P> &v, Args... args);
	template <unsigned int X, typename... Args>
	void initialize(const T &v, Args... args);
	
public:
	// Catch-all variadic constructor
	// The following constructor calls are supported:
	// Matrix(): Creates an empty matrix
	//  All elements are 0.
	// Matrix(const T &d): Identity or scaling matrix
	//  The diagonal elements are all d, everything else is 0.
	// Matrix<ForwardIterator>(std::pair<ForwardIterator, ForwardIterator> iterator): Populate the matrix from an iterator pair
	//  Elements are filled up in order of occurrence, column by column, row by row. Transpose if you want the same layout as in a const array.
	//  std::pair is used to avoid call ambiguity problems. Call std::make_pair(begin, end) to constuct your iterator pair.
	// Matrix<TypeName, Height, Width>(const Matrix<TypeName, Height, Width> &m): Create a plain, cast or resized copy of another matrix
	//  Elements are filled up column by column, padded up with 0s. The diagonal is padded with 1s.
	// Matrix<AnyScalarOrVector...>(const AnyScalarOrVector &arg0, ...): Populate the matrix from a list of vectors and/or scalars
	//  Elements are filled up in order of occurrence, column by column, row by row.
	// Matrix(const T *m): Populate the matrix from a constant array
	//  Elements are filled up in order of occurrence, column by column, row by row. Transpose if you want the same layout as in a const array.
	//  Note that this overload does not support range checking. Make sure that the input array has the correct size or use the ForwardIterator overload.
	template <typename... Args>
	Matrix(Args... args);
	
	// Rotation matrices
	// Angle must be specified in radians
	// Only implemented for mat2 (2D linear transform) and Matrix3f (2D affine transform)
	static Matrix<T, M, N> rotation(const T &angle) __attribute__((deprecated));
	// Only implemented for mat3 and mat4 (3D linear/affine transform)
	static Matrix<T, M, N> rotation(const T &angle, const T &a, const T &b, const T &c) __attribute__((deprecated));
	// Only implemented for mat3 (3D linear transform)
	static Matrix<T, M, N> rotation(const T &angle, const Vector<T, M> &v) __attribute__((deprecated));
	// Only implemented for mat4 (3D affine transform)
	static Matrix<T, M, N> rotation(const T &angle, const Vector<T, M - 1> &v) __attribute__((deprecated));
	// Scaling matrix
	// Only implemented for M = N
	static Matrix<T, M, N> scaling(const Vector<T, M> &v) __attribute__((deprecated));
	// Only implemented for mat3 and mat4 (3D linear/affine transform)
	static Matrix<T, M, N> scaling(const T &a, const T &b, const T &c) __attribute__((deprecated));
	// Affine translation matrix
	// Only implemented for M = N
	static Matrix<T, M, N> translation(const Vector<T, M - 1> &v) __attribute__((deprecated));
	// Only implemented for mat3 and mat4 (affine transform)
	static Matrix<T, M, N> translation(const T &a, const T &b, const T &c = 0) __attribute__((deprecated));
	// Perspective projection matrix, specification through field-of-view angle, aspect ratio and near/far planes
	// Angle must be specified in radians
	// Only implemented for mat4
	static Matrix<T, M, N> perspective(const T &fovy, const T &aspect, const T &nearz, const T &farz) __attribute__((deprecated));
	// Perspective projection matrix, specification through view frustum edges (left, right, bottom, top, near, far)
	// Only implemented for mat4
	static Matrix<T, M, N> frustum(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f) __attribute__((deprecated));
	// Orthographic projection matrix, specification through view cube edges (left, right, bottom, top, near, far)
	// Only implemented for mat4
	static Matrix<T, M, N> ortho(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f) __attribute__((deprecated));

	// Get internal pointer, subject to memory alignment (padding between column vectors)
	const T *internal() const;
	// Column vector access operator
	// Use double index operators (i.e. [j][i]) to access a single element
	Vector<T, M> &operator [](unsigned int j);
	// Const column vector access operator
	// Use double index operators (i.e. [j][i]) to access a single element
	const Vector<T, M> &operator [](unsigned int j) const;
	// Assignment operator
	Matrix<T, M, N> &operator =(const Matrix<T, M, N> &other);
	// Component-wise sum
	Matrix<T, M, N> &operator +=(const Matrix<T, M, N> &other);
	// Component-wise difference
	Matrix<T, M, N> &operator -=(const Matrix<T, M, N> &other);
	// Matrix-scalar product
	Matrix<T, M, N> &operator *=(const T &x);
	// Matrix-scalar division
	Matrix<T, M, N> &operator /=(const T &x);

	// Transpose
	Matrix<T, N, M> t() const __attribute__((deprecated));

	// Matrix determinant through LU decomposition
	T det() const __attribute__((deprecated));
	// Matrix inverse through LU decomposition
	Matrix<T, M, N> inv() const __attribute__((deprecated));
};

// Component-wise matrix sum
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b);

// Component-wise matrix difference
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b);

// Multiplication of every matrix element with a factor
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator *(const Matrix<T, M, N> &a, const T &x);

// Division of every matrix element by a divisor
template <class T, unsigned int M, unsigned int N, class U>
inline Matrix<T, M, N> operator /(const Matrix<T, M, N> &a, const U &x);

// Matrix-vector product
template <class T, unsigned int M, unsigned int N>
inline Vector<T, N> operator *(const Matrix<T, M, N> &a, const Vector<T, N> &b);

// Matrix product
template <class T, unsigned int M, unsigned int N, unsigned int P>
inline Matrix<T, M, P> operator *(const Matrix<T, M, N> &a, const Matrix<T, N, P> &b);

// Matrix comparison operator (exact comparison)
template <class T, unsigned int M, unsigned int N>
inline bool operator ==(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b);

// Component-wise matrix product
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> matrixCompMult(const Matrix<T, M, N> &x, const Matrix<T, M, N> &y);

// Cartesian vector product
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> outerProduct(const Vector<T, N> &c, const Vector<T, M> &r);

// Matrix transpose (columns and rows are swapped)
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, N, M> transpose(const Matrix<T, M, N> &m);

// Matrix determinant
template <class T, unsigned int M, unsigned int N>
inline T determinant(const Matrix<T, M, N> &m);

// Matrix inverse (through LU decomposition)
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> inverse(const Matrix<T, M, N> &m);

// Output stream operator
template <class T, unsigned int M, unsigned int N>
inline std::ostream &operator <<(std::ostream &o, const Matrix<T, M, N> &m);

// Generate a rotation matrix
// Angle must be specified in radians.
// 2D linear transform
template <class T>
inline Matrix<T, 2, 2> rotationMatrix(const T &angle);
// 3D linear transform
template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const T &x, const T &y, const T &z);
template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const Vector<T, 3> &v);

// Generate a scaling matrix
// 2D linear transform
template <class T>
inline Matrix<T, 2, 2> scalingMatrix(const T &x, const T &y);
// 3D linear transform
template <class T>
inline Matrix<T, 3, 3> scalingMatrix(const T &x, const T &y, const T &z);
// 4D linear transform
template <class T>
inline Matrix<T, 4, 4> scalingMatrix(const T &x, const T &y, const T &z, const T &w);
// M-dimensional linear transform
template <class T, unsigned int M>
inline Matrix<T, M, M> scalingMatrix(const Vector<T, M> &v);

// Generate a translation matrix
// 2D affine transform
template <class T>
inline Matrix<T, 3, 3> translationMatrix(const Vector<T, 2> &v);
template <class T>
inline Matrix<T, 3, 3> translationMatrix(const T &x, const T &y);
// 3D affine transform
template <class T>
inline Matrix<T, 4, 4> translationMatrix(const Vector<T, 3> &v);
template <class T>
inline Matrix<T, 4, 4> translationMatrix(const T &x, const T &y, const T &z);

// Generate a perspective projection matrix
// Specification through field-of-view angle, aspect ratio and near/far planes.
// Angle must be specified in radians.
template <class T>
inline Matrix<T, 4, 4> perspectiveMatrix(const T &fovy, const T &aspect, const T &nearz, const T &farz);
// Specification through view frustum edges (left, right, bottom, top, near, far)
template <class T>
inline Matrix<T, 4, 4> frustumMatrix(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f);

// Generate an orthographic projection matrix
// Specification through view cube edges (left, right, bottom, top, near, far).
template <class T>
inline Matrix<T, 4, 4> orthoMatrix(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f);

// Predefined types
typedef Matrix<float, 2, 2> mat2;
typedef Matrix<float, 3, 3> mat3;
typedef Matrix<float, 4, 4> mat4;
typedef Matrix<float, 2, 2> mat2x2;
typedef Matrix<float, 3, 2> mat2x3;
typedef Matrix<float, 4, 2> mat2x4;
typedef Matrix<float, 2, 3> mat3x2;
typedef Matrix<float, 3, 3> mat3x3;
typedef Matrix<float, 4, 3> mat3x4;
typedef Matrix<float, 2, 4> mat4x2;
typedef Matrix<float, 3, 4> mat4x3;
typedef Matrix<float, 4, 4> mat4x4;

// Implementation

#warning TODO Remove deprecation pragma when everyting is moved to the function API
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

template <class T, unsigned int M, unsigned int N>
T &Matrix<T, M, N>::element(unsigned int index) {
	return _a[index / M][index % M];
}

template <class T, unsigned int M, unsigned int N>
const T &Matrix<T, M, N>::element(unsigned int index) const {
	return _a[index / M][index % M];
}

template <class T, unsigned int M, unsigned int N>
template <unsigned int X>
inline void Matrix<T, M, N>::initialize() { }

template <class T, unsigned int M, unsigned int N>
template <unsigned int X, typename... Args>
inline void Matrix<T, M, N>::initialize(const T &v, Args... args) {
	// The condition is < and not <= because there is a dedicated terminal initializer for a single scalar
	static_assert(X + 1 < M * N, "Too many initializers");
	element(X) = v;
	initialize<X + 1>(args...);
}

template <class T, unsigned int M, unsigned int N>
template <unsigned int X, typename... Args, unsigned int P>
inline void Matrix<T, M, N>::initialize(const Vector<T, P> &v, Args... args) {
	static_assert(X + P <= M * N, "Too many initializers");
	for (unsigned int p = 0; p < P; p++) {
		element(X + p) = v[p];
	}
	initialize<X + P>(args...);
}

template <class T, unsigned int M, unsigned int N>
template <unsigned int X>
inline void Matrix<T, M, N>::initialize(const T &v) {
	static_assert(X + 1 <= M * N, "Too many initializers");
	if (X == 0) {
		for (unsigned int i = 0; i < M; i++) {
			for (unsigned int j = 0; j < N; j++) {
				if (i == j) {
					// Diagonal element, fill with value
					(*this)[j][i] = v;
				} else {
					// Anything else, fill with 0
					(*this)[j][i] = static_cast<const T>(0);
				}
			}
		}
	} else {
		element(X) = v;
	}
}

template <class T, unsigned int M, unsigned int N>
template <unsigned int X, typename ForwardIterator>
inline void Matrix<T, M, N>::initialize(std::pair<ForwardIterator, ForwardIterator> iterator) {
	static_assert(X == 0, "No other arguments are allowed when initializing a matrix from an iterator");
	unsigned int index = X;
	for (ForwardIterator it = iterator.first; it != iterator.second; it++, index++) {
		element(index) = *it;
	}
}

template <class T, unsigned int M, unsigned int N>
template <unsigned int X>
inline void Matrix<T, M, N>::initialize(const T *v) {
	static_assert(X == 0, "No other arguments are allowed when initializing a matrix from an array");
	for (unsigned int p = 0; p < M * N; p++) {
		element(X + p) = v[p];
	}
}

template <class T, unsigned int M, unsigned int N>
template <unsigned int X, class U, unsigned int P, unsigned int R>
inline void Matrix<T, M, N>::initialize(const Matrix<U, P, R> &m) {
	static_assert(X == 0, "No other arguments are allowed when copy-constructing a matrix");
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			if (i < P && j < R) {
				// Component available, copy
				(*this)[j][i] = static_cast<const T>(m[j][i]);
			} else {
				// Component not available, fill
				if (i == j) {
					// Diagonal element, fill with identity
					(*this)[j][i] = static_cast<const T>(1);
				} else {
					// Anything else, fill with 0
					(*this)[j][i] = static_cast<const T>(0);
				}
			}
		}
	}
}

template <class T, unsigned int M, unsigned int N>
template <typename... Args>
inline Matrix<T, M, N>::Matrix(Args... args) {
	initialize<0>(args...);
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator =(const Matrix<T, M, N> &other) {
	for (unsigned int j = 0; j < N; j++) {
		_a[j] = other[j];
	}
	return *this;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator +=(const Matrix<T, M, N> &other) {
	for (unsigned int j = 0; j < N; j++) {
		_a[j] += other[j];
	}
	return *this;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b) {
	Matrix<T, M, N> ret = a;
	ret += b;
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator -=(const Matrix<T, M, N> &other) {
	for (unsigned int j = 0; j < N; j++) {
		_a[j] -= other[j];
	}
	return *this;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b) {
	Matrix<T, M, N> ret = a;
	ret -= b;
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator *=(const T &x) {
	for (unsigned int i = 0; i < M * N; i++) {
		_a[i] *= x;
	}
	return *this;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator *(const Matrix<T, M, N> &a, const T &x) {
	Matrix<T, M, N> ret = a;
	ret *= x;
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator /=(const T &x) {
	*this *= 1 / x;
	return *this;
}

template <class T, unsigned int M, unsigned int N, class U>
inline Matrix<T, M, N> operator /(const Matrix<T, M, N> &a, const U &x) {
	Matrix<T, M, N> ret = a;
	ret /= x;
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline std::ostream &operator <<(std::ostream &o, const Matrix<T, M, N> &m) {
	o << "{";
	for (unsigned int i = 0; i < M; i++) {
		o << "\n\t";
		for (unsigned int j = 0; j < N; j++) {
			o << m[j][i] << "\t";
		}
	}
	o << "\n}";
	return o;
}

template <class T, unsigned int M, unsigned int N>
inline bool operator ==(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b) {
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			if (a[j][i] != b[j][i]) {
				return false;
			}
		}
	}
	return true;
}

template <class T, unsigned int M, unsigned int N, unsigned int P>
inline Matrix<T, M, P> operator *(const Matrix<T, M, N> &a, const Matrix<T, N, P> &b) {
	//      4 x 3               3 x 4              4 x 4
	// | a00 a01 a02 |   | b00 b01 b02 b03 |   | a00*b00+a01*b10+a02*b20 a00*b01+a01*b11+a02*b21 ... |
	// | a10 a11 a12 | * | b10 b11 b12 b33 | = | a10*b00+a11*b10+a12*b20 ...                         |
	// | a20 a21 a22 |   | b20 b21 b22 b33 |   | ...                                                 |
	// | a30 a31 a32 |                         | ...                                                 |
	Matrix<T, M, P> ret;
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < P; j++) {
			ret[j][i] = 0;
			for (unsigned int k = 0; k < N; k++) {
				ret[j][i] += a[k][i] * b[j][k];
			}
		}
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Vector<T, N> operator *(const Matrix<T, M, N> &a, const Vector<T, N> &b) {
	Vector<T, M> ret;
	for (unsigned int i = 0; i < M; i++) {
		ret[i] = 0;
		for (unsigned int k = 0; k < N; k++) {
			ret[i] += a[k][i] * b[k];
		}
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline const T *Matrix<T, M, N>::internal() const {
	return &(*this)[0][0];
}

template <class T, unsigned int M, unsigned int N>
inline Vector<T, M> &Matrix<T, M, N>::operator [](unsigned int j) {
	if (j < N) {
		return _a[j];
	}
	throw DimensionOutOfRangeException("Matrix column index out of range", N, j);
}

template <class T, unsigned int M, unsigned int N>
inline const Vector<T, M> &Matrix<T, M, N>::operator [](unsigned int j) const {
	if (j < N) {
		return _a[j];
	}
	throw DimensionOutOfRangeException("Matrix column index out of range", N, j);
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, N, M> Matrix<T, M, N>::t() const {
	Matrix<T, N, M> ret;
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			ret[i][j] = (*this)[j][i];
		}
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline std::tuple<Matrix<T, M, N>, Matrix<T, M, N>, Matrix<T, M, N>, unsigned int> Matrix<T, M, N>::lu() const {
	static_assert(M == N, "Matrix is not square");
	// Prepare result (lower and upper triangular matrices, permutation)
	Matrix<T, M, N> l(1);
	Matrix<T, M, N> u = *this;
	Matrix<T, M, N> p(1);
	// Reset the swap counter
	unsigned int swaps = 0;
	// Loop through all rows except the last (which will simply become 0 0 .. 0 1 in L)
	for (unsigned int n = 0; n < M - 1; n++) {
		// Prepare the intermediate lower matrix (one column filled)
		Matrix<T, M, N> ln(1);
		// Check if a swap is needed
		unsigned int s = M;
		for (unsigned int i = n; i < M && s == M; i++) {
			if (u[i][n] != 0) {
				s = i;
			}
		}
		if (s == M) {
			// All coefficients of column n are 0
			throw NonInvertibleMatrixException("Can't find solution for matrix inverse");
		} else if (s != n) {
			// (n, n) is 0, swap row n with row s (which has a non-zero coefficient)
			// Also swap the corresponding rows in the permutation matrix
			for (unsigned int j = 0; j < N; j++) {
				std::swap(u[j][n], u[j][s]);
				std::swap(p[j][n], p[j][s]);
			}
			// Increment the swap counter
			swaps++;
		}
		// Calculate lower matrix coefficients for this column
		for (unsigned int i = n + 1; i < M; i++) {
			T v = u[n][i] / u[n][n];
			// And assign them to the appropriate Ln
			ln[n][i] = -v;
			// and L fields
			l[n][i] = v;
		}
		// Apply Ln to U, this cancels out all coefficients under U(n,n)
		u = ln * u;
	}
	return std::make_tuple(l, u, p, swaps);
}


template <class T, unsigned int M, unsigned int N>
T Matrix<T, M, N>::det() const {
	static_assert(M == N && M > 0, "Matrix is not square");
	// Decompose into triangular parts
	try {
		std::tuple<Matrix<T, M, N>, Matrix<T, M, N>, Matrix<T, M, N>, unsigned int> lups = lu();
		// det(A) = det(P^-1) * det(L) * det(U) = (-1)^S * (l11 * l22 * ...) * (u11 * u22 * ...)
		// Since the diagonal of L is all 1s, its determinant is also 1, so
		// det(A) = (-1)^S * (u11 * u22 * ...)
		// Calculate d according to (-1)^S here, S = number of exchanges in P^-1
		T d;
		if (std::get<SWAPS>(lups) & 1) {
			d = static_cast<T>(-1);
		} else {
			d = static_cast<T>(1);
		}
		for (unsigned int i = 0; i < M; i++) {
			d *= std::get<UPPER>(lups)[i][i];
		}
		return d;
	} catch (NonInvertibleMatrixException &ex) {
		return static_cast<T>(0);
	}
}

template <class T, unsigned int M, unsigned int N>
Matrix<T, M, N> Matrix<T, M, N>::inv() const {
	static_assert(M == N && M > 0, "Matrix is not square");
	// Decompose into triangular parts
	std::tuple<Matrix<T, M, N>, Matrix<T, M, N>, Matrix<T, M, N>, unsigned int> lups = lu();
	// Start with AX = I, where X is the inverse of A
	// A = P^-1LU, so P^-1LUX = I, and thus LUX = P (with P = I if no row swapping was needed)
	// Substitute UX = Y, yielding LY = P
	// Calculate Y through forward substitution of L
	Matrix<T, M, N> y;
	// Process each column independently
	for (unsigned int j = 0; j < N; j++) {
		// Substitute each row from the previous rows
		for (unsigned int i = 0; i < M; i++) {
			// Fetch the starting value from P
			T diff = std::get<PERMUTATION>(lups)[j][i];
			// Subtract the product of the values from the previous rows (from the top) times the corresponding values in L
			for (unsigned int k = 0; k < i; k++) {
				diff -= std::get<LOWER>(lups)[k][i] * y[j][k];
			}
			// No division necessary, the diagonal is always 1
			y[j][i] = diff;
		}
	}
	// Calculate X from UX = Y through backward substitution of U
	Matrix<T, M, N> x;
	// Process each column independently
	for (unsigned int j = 0; j < N; j++) {
		// Substitute each row from the previous rows, work backwards
		for (unsigned int i = M; i-- > 0;) {
			// Fetch the starting value from Y
			T diff = y[j][i];
			// Subtract the product of the values from the previous rows (from the bottom) times the corresponding values in U
			for (unsigned int k = i + 1; k < M; k++) {
				diff -= std::get<UPPER>(lups)[k][i] * x[j][k];
			}
			// Divide by the diagonal coefficient
			x[j][i] = diff / std::get<UPPER>(lups)[i][i];
		}
	}
	// X is the inverse of A
	return x;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> matrixCompMult(const Matrix<T, M, N> &x, const Matrix<T, M, N> &y) {
	Matrix<T, M, N> ret;
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			ret[j][i] = x[j][i] * y[j][i];
		}
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> outerProduct(const Vector<T, M> &c, const Vector<T, N> &r) {
	return Matrix<T, M, 1>(c) * Matrix<T, 1, N>(r);
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, N, M> transpose(const Matrix<T, M, N> &m) {
	return m.t();
}

template <class T, unsigned int M, unsigned int N>
inline T determinant(const Matrix<T, M, N> &m){
	return m.det();
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> inverse(const Matrix<T, M, N> &m) {
	return m.inv();
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> Matrix<T, M, N>::scaling(const Vector<T, M> &v) {
	static_assert(M == N && M > 0, "Matrix is not square");
	Matrix<T, M, M> ret(1);
	for (unsigned int i = 0; i < M; i++) {
		ret[i][i] = v[i];
	}
	return ret;
}

template <>
inline Matrix<float, 3, 3> Matrix<float, 3, 3>::scaling(const float &a, const float &b, const float &c) {
	return Matrix<float, 3, 3>::scaling(Vector<float, 3>(a, b, c));
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::scaling(const float &a, const float &b, const float &c) {
	return Matrix<float, 4, 4>(Matrix<float, 3, 3>::scaling(Vector<float, 3>(a, b, c)));
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> Matrix<T, M, N>::translation(const Vector<T, M - 1> &v) {
	static_assert(M == N && M > 0, "Matrix is not square");
	Matrix<T, M, M> ret(1);
	for (unsigned int i = 0; i < M - 1; i++) {
		ret[N - 1][i] = v[i];
	}
	// As long as the conversion constructor is explicit, this will fail if the matrix isn't square
	return ret;
}

template <>
inline Matrix<float, 3, 3> Matrix<float, 3, 3>::translation(const float &a, const float &b, const float &c) {
	return Matrix<float, 3, 3>::translation(Vector<float, 2>(a, b));
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::translation(const float &a, const float &b, const float &c) {
	return Matrix<float, 4, 4>::translation(Vector<float, 3>(a, b, c));
}

template <>
inline Matrix<float, 2, 2> Matrix<float, 2, 2>::rotation(const float &angle) {
	float s = sin(angle);
	float c = cos(angle);
	// note column-major order
	float data[] = {
		c, s,
		-s, c,
	};
	return Matrix<float, 2, 2>(data);
}

template <>
inline Matrix<float, 3, 3> Matrix<float, 3, 3>::rotation(const float &angle, const Vector<float, 3> &v) {
	Vector<float, 3> u = normalize(v);
	float s = sin(angle);
	float c = cos(angle);
	// note column-major order
	float data[] = {
		u[0] * u[0] + (1 - u[0] * u[0]) * c, u[0] * u[1] * (1 - c) + u[2] * s, u[0] * u[2] * (1 - c) - u[1] * s,
		u[0] * u[1] * (1 - c) - u[2] * s, u[1] * u[1] + (1 - u[1] * u[1]) * c, u[1] * u[2] * (1 - c) + u[0] * s,
		u[0] * u[2] * (1 - c) + u[1] * s, u[1] * u[2] * (1 - c) - u[0] * s, u[2] * u[2] + (1 - u[2] * u[2]) * c,
	};
	return Matrix<float, 3, 3>(data);
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::rotation(const float &angle, const Vector<float, 3> &v) {
	return Matrix<float, 4, 4>(Matrix<float, 3, 3>::rotation(angle, v));
}

template <>
inline Matrix<float, 3, 3> Matrix<float, 3, 3>::rotation(const float &angle, const float &a, const float &b, const float &c) {
	return Matrix<float, 3, 3>::rotation(angle, Vector<float, 3>(a, b, c));
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::rotation(const float &angle, const float &a, const float &b, const float &c) {
	return Matrix<float, 4, 4>::rotation(angle, Vector<float, 3>(a, b, c));
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::frustum(const float &l, const float &r, const float &b, const float &t, const float &n, const float &f) {
	float data[] = {
		2.0f * n / (r - l), 0, 0, 0,
		0, 2.0f * n / (t - b), 0, 0,
		(r + l) / (r - l), (t + b) / (t - b), -((f + n) / (f - n)), -1.0f,
		0, 0, -2.0f * f * n / (f - n), 0,
	};
	return Matrix<float, 4, 4>(data);
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::perspective(const float &fovy, const float &aspect, const float &nearz, const float &farz) {
	float f = 1.0f / tan(fovy / 2.0f);
	// note column-major order
	float data[] = {
		f / aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (farz + nearz) / (nearz - farz), -1.0f,
		0, 0, 2.0f * farz * nearz / (nearz - farz), 0
	};
	return Matrix<float, 4, 4>(data);
}

template <>
inline Matrix<float, 4, 4> Matrix<float, 4, 4>::ortho(const float &l, const float &r, const float &b, const float &t, const float &n, const float &f) {
	float data[] = {
		2.0f / (r - l), 0, 0, 0,
		0, 2.0f / (t - b), 0, 0,
		0, 0, -2.0f / (f - n), 0,
		-((r + l) / (r - l)), -((t + b) / (t - b)), -((f + n) / (f - n)), 1,
	};
	return Matrix<float, 4, 4>(data);
}

template <class T>
inline Matrix<T, 2, 2> rotationMatrix(const T &angle) {
	return Matrix<T, 2, 2>::rotation(angle);
}

template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle) {
	return Matrix<T, 3, 3>::rotation(angle);
}

template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const T &x, const T &y, const T &z){
	return Matrix<T, 3, 3>::rotation(angle, x, y, z);
}

template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const Vector<T, 3> &v) {
	return Matrix<T, 3, 3>::rotation(angle, v);
}

template <class T>
inline Matrix<T, 2, 2> scalingMatrix(const T &x, const T &y) {
	return Matrix<T, 2, 2>::scaling(x, y);
}

template <class T>
inline Matrix<T, 3, 3> scalingMatrix(const T &x, const T &y, const T &z) {
	return Matrix<T, 3, 3>::scaling(x, y, z);
}

template <class T>
inline Matrix<T, 4, 4> scalingMatrix(const T &x, const T &y, const T &z, const T &w) {
	return Matrix<T, 2, 2>::scaling(x, y, z, w);
}

template <class T, unsigned int M>
inline Matrix<T, M, M> scalingMatrix(const Vector<T, M> &v) {
	return Matrix<T, M, M>::scaling(v);
}

template <class T>
inline Matrix<T, 3, 3> translationMatrix(const Vector<T, 2> &v) {
	return Matrix<T, 3, 3>::translation(v);
}

template <class T>
inline Matrix<T, 3, 3> translationMatrix(const T &x, const T &y) {
	return Matrix<T, 3, 3>::translation(x, y);
}

template <class T>
inline Matrix<T, 4, 4> translationMatrix(const Vector<T, 3> &v) {
	return Matrix<T, 4, 4>::translation(v);
}

template <class T>
inline Matrix<T, 4, 4> translationMatrix(const T &x, const T &y, const T &z) {
	return Matrix<T, 4, 4>::translation(x, y, z);
}

template <class T>
inline Matrix<T, 4, 4> perspectiveMatrix(const T &fovy, const T &aspect, const T &nearz, const T &farz) {
	return Matrix<T, 4, 4>::perspective(fovy, aspect, nearz, farz);
}

template <class T>
inline Matrix<T, 4, 4> frustumMatrix(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f) {
	return Matrix<T, 4, 4>::frustum(l, r, b, t, n, f);
}

template <class T>
inline Matrix<T, 4, 4> orthoMatrix(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f) {
	return Matrix<T, 4, 4>::ortho(l, r, b, t, n, f);
}

#warning TODO Remove deprecation pragma when everyting is moved to the function API
#pragma GCC diagnostic pop

}

#endif //GLAM_MATRIX_H

