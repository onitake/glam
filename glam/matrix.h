/*
 * GLAM - GLSL Linear Algebra Math Library
 * 
 * Copyright (c) 2012-2014, Gregor Riepl <onitake@gmail.com>
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
#include <glam/vector.h>
#include <glam/exception.h>
#include <glam/config.h>

namespace glam {

// M=height N=width i=row j=column
template <class T, unsigned int M, unsigned int N = M>
class Matrix {
private:
	// The matrix is sorted in column-major order, as in OpenGL:
	// | a00 a01 a02 | => [ a00 a10 ] [ a01 a11 ] [ a02 a12 ] (index = M * j + i)
	// | a10 a11 a12 |
	Vector<T, M> _a[N];

	// Single element access by array index (inefficient implementation, use sparingly)
	T &element(unsigned int index);
	inline const T &element(unsigned int index) const;

	// Supporting initializers to avoid trouble from delegate constructors
	template <unsigned int X>
	inline void initialize();
	template <unsigned int X>
	inline void initialize(const T &v);
	template <unsigned int X, typename ForwardIterator>
	inline void initialize(std::pair<ForwardIterator, ForwardIterator> iterator);
	template <unsigned int X>
	inline void initialize(const T *v);
	template <unsigned int X, class U, unsigned int P, unsigned int R>
	inline void initialize(const Matrix<U, P, R> &m);
	template <unsigned int X, typename... Args, unsigned int P>
	inline void initialize(const Vector<T, P> &v, Args... args);
	template <unsigned int X, typename... Args>
	inline void initialize(const T &v, Args... args);
	
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
	inline Matrix(Args... args);
	
	// Get internal pointer, subject to memory alignment (padding between column vectors)
	inline const T *internal() const;
	// Column vector access operator
	// Use double index operators (i.e. [j][i]) to access a single element
	inline Vector<T, M> &operator [](unsigned int j);
	// Const column vector access operator
	// Use double index operators (i.e. [j][i]) to access a single element
	inline const Vector<T, M> &operator [](unsigned int j) const;
	// Assignment operator
	inline Matrix<T, M, N> &operator =(const Matrix<T, M, N> &other);
	// Component-wise sum
	inline Matrix<T, M, N> &operator +=(const Matrix<T, M, N> &other);
	// Matrix-scalar sum
	inline Matrix<T, M, N> &operator +=(const T &x);
	// Component-wise difference
	inline 	Matrix<T, M, N> &operator -=(const Matrix<T, M, N> &other);
	// Matrix-scalar difference
	inline Matrix<T, M, N> &operator -=(const T &x);
	// Matrix division
	inline Matrix<T, M, N> &operator /=(const Matrix<T, M, N> &other);
	// Matrix-scalar division
	inline Matrix<T, M, N> &operator /=(const T &x);
	// Matrix-scalar product
	inline Matrix<T, M, N> &operator *=(const T &x);
};

// LU(P) decomposition of M
// lower is a lower triangular matrix (L)
// upper is an upper triangular matrix (U)
// permutation is a permutation matrix (P)
// swaps is the number of swaps performed to obtain permutation, with respect to identity (S)
// The following equation must hold: PM = LU
template <class T, unsigned int M, unsigned int N>
struct LuDecomposition {
	Matrix<T, M, N> lower;
	Matrix<T, M, N> upper;
	Matrix<T, M, N> permutation;
	unsigned int swaps;
	inline LuDecomposition(const Matrix<T, M, N> &m);
};

// Component-wise matrix sum
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b);
// Component-wise matrix-scalar sum
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const Matrix<T, M, N> &a, const T &b);
// Component-wise scalar-matrix sum
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const T &a, const Matrix<T, M, N> &b);

// Component-wise matrix difference
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b);
// Component-wise matrix-scalar difference
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const Matrix<T, M, N> &a, const T &b);
// Component-wise scalar-matrix difference
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const T &a, const Matrix<T, M, N> &b);

// Component-wise matrix division
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator /(const Matrix<T, M, N> &a, const Matrix<T, M, N> &b);
// Division of every matrix component by a scalar
template <class T, unsigned int M, unsigned int N, class U>
inline Matrix<T, M, N> operator /(const Matrix<T, M, N> &a, const U &x);
// Division of a scalar with every matrix component
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator /(const T &x, const Matrix<T, M, N> &a);

// Multiplication of every matrix component with a scalar
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator *(const Matrix<T, M, N> &a, const T &x);
// Multiplication of a scalar with every matrix component
template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator *(const T &x, const Matrix<T, M, N> &a);

// Matrix-vector product
template <class T, unsigned int M, unsigned int N>
inline Vector<T, M> operator *(const Matrix<T, M, N> &a, const Vector<T, N> &b);
// Vector-matrix product
template <class T, unsigned int M, unsigned int N>
inline Vector<T, M> operator *(const Vector<T, N> &a, const Matrix<T, M, N> &b);

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
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const Vector<T, 3> &v);
template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const T &x, const T &y, const T &z);

// Generate a scaling matrix
// Linear transform
template <class T, unsigned int M>
inline Matrix<T, M, M> scalingMatrix(const Vector<T, M> &v);
// 2D linear transform
template <class T>
inline Matrix<T, 2, 2> scalingMatrix(const T &x, const T &y);
// 3D linear transform
template <class T>
inline Matrix<T, 3, 3> scalingMatrix(const T &x, const T &y, const T &z);

// Generate a translation matrix
// Affine transform
template <class T, unsigned int M>
inline Matrix<T, M + 1, M + 1> translationMatrix(const Vector<T, M> &v);
// 2D affine transform
template <class T>
inline Matrix<T, 3, 3> translationMatrix(const T &x, const T &y);
// 3D affine transform
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

template <class T, unsigned int M, unsigned int N>
inline T &Matrix<T, M, N>::element(unsigned int index) {
	return (*this)[index / M][index % M];
}

template <class T, unsigned int M, unsigned int N>
inline const T &Matrix<T, M, N>::element(unsigned int index) const {
	return (*this)[index / M][index % M];
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
	for (unsigned int j = 0; j < N; j++) {
		if (a[j] != b[j]) {
			return false;
		}
	}
	/*for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			if (a[j][i] != b[j][i]) {
				return false;
			}
		}
	}*/
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
	Matrix<T, N, M> t = transpose(a);
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < P; j++) {
			ret[j][i] = dot(t[i], b[j]);
		}
	}
	/*for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < P; j++) {
			ret[j][i] = 0;
			for (unsigned int k = 0; k < N; k++) {
				ret[j][i] += a[k][i] * b[j][k];
			}
		}
	}*/
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Vector<T, M> operator *(const Matrix<T, M, N> &a, const Vector<T, N> &b) {
	// Transpose the matrix
	Vector<T, M> ret;
	Matrix<T, N, M> t = transpose(a);
	for (unsigned int i = 0; i < M; i++) {
		ret[i] = dot(t[i], b);
	}
	/*for (unsigned int i = 0; i < M; i++) {
		ret[i] = 0;
		for (unsigned int k = 0; k < N; k++) {
			ret[i] += a[k][i] * b[k];
		}
	}*/
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline const T *Matrix<T, M, N>::internal() const {
	return &(*this)[0][0];
}

template <class T, unsigned int M, unsigned int N>
inline Vector<T, M> &Matrix<T, M, N>::operator [](unsigned int j) {
#ifdef GLAM_RANGE_CHECKS
	if (j < N) {
		return _a[j];
	}
	throw DimensionOutOfRangeException("Matrix column index out of range", N, j);
#else
	return _a[j];
#endif
}

template <class T, unsigned int M, unsigned int N>
inline const Vector<T, M> &Matrix<T, M, N>::operator [](unsigned int j) const {
#ifdef GLAM_RANGE_CHECKS
	if (j < N) {
		return _a[j];
	}
	throw DimensionOutOfRangeException("Matrix column index out of range", N, j);
#else
	return _a[j];
#endif
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
inline LuDecomposition<T, M, N>::LuDecomposition(const Matrix<T, M, N> &m) : lower(1), upper(m), permutation(1), swaps(0) {
	static_assert(M == N, "Matrix is not square");
	// Initializer list: prepare result (lower and upper triangular matrices, permutation, reset the swap counter)
	// Loop through all rows except the last (which will simply become 0 0 .. 0 1 in L)
	for (unsigned int n = 0; n < M - 1; n++) {
		// Prepare the intermediate lower matrix (one column filled)
		Matrix<T, M, N> ln(1);
		// Check if a swap is needed
		unsigned int s = M;
		for (unsigned int i = n; i < M && s == M; i++) {
			if (upper[i][n] != 0) {
				s = i;
			}
		}
		if (s == M) {
			// All coefficients of column n are 0
			throw NonInvertibleMatrixException("Singular matrix, can't find solution for inverse");
		} else if (s != n) {
			// (n, n) is 0, swap row n with row s (which has a non-zero coefficient)
			// Also swap the corresponding rows in the permutation matrix
			for (unsigned int j = 0; j < N; j++) {
				std::swap(upper[j][n], upper[j][s]);
				std::swap(permutation[j][n], permutation[j][s]);
			}
			// Increment the swap counter
			swaps++;
		}
		// Calculate lower matrix coefficients for this column
		for (unsigned int i = n + 1; i < M; i++) {
			T v = upper[n][i] / upper[n][n];
			// And assign them to the appropriate Ln
			ln[n][i] = -v;
			// and L fields
			lower[n][i] = v;
		}
		// Apply Ln to U, this cancels out all coefficients under U(n,n)
		upper = ln * upper;
	}
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, N, M> transpose(const Matrix<T, M, N> &m) {
	Matrix<T, N, M> ret;
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < N; j++) {
			ret[i][j] = m[j][i];
		}
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline T determinant(const Matrix<T, M, N> &m) {
	static_assert(M == N && M > 0, "Matrix is not square");
	// Decompose into triangular parts
	try {
		LuDecomposition<T, M, N> lups(m);
		// det(A) = det(P^-1) * det(L) * det(U) = (-1)^S * (l11 * l22 * ...) * (u11 * u22 * ...)
		// Since the diagonal of L is all 1s, its determinant is also 1, so
		// det(A) = (-1)^S * (u11 * u22 * ...)
		// Calculate d according to (-1)^S here, S = number of exchanges in P^-1
		T d;
		if (lups.swaps & 1) {
			d = static_cast<T>(-1);
		} else {
			d = static_cast<T>(1);
		}
		for (unsigned int i = 0; i < M; i++) {
			d *= lups.upper[i][i];
		}
		return d;
	} catch (NonInvertibleMatrixException &ex) {
		return T(0);
	}
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> inverse(const Matrix<T, M, N> &m) {
	static_assert(M == N && M > 0, "Matrix is not square");
	// Decompose into triangular parts
	LuDecomposition<T, M, N> lups(m);
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
			T diff = lups.permutation[j][i];
			// Subtract the product of the values from the previous rows (from the top) times the corresponding values in L
			for (unsigned int k = 0; k < i; k++) {
				diff -= lups.lower[k][i] * y[j][k];
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
				diff -= lups.upper[k][i] * x[j][k];
			}
			// Divide by the diagonal coefficient
			x[j][i] = diff / lups.upper[i][i];
		}
	}
	// X is the inverse of A
	return x;
}

template <class T>
inline Matrix<T, 2, 2> rotationMatrix(const T &angle) {
	T s = sin(angle);
	T c = cos(angle);
	// note column-major order
	T data[] = {
		c, s,
		-s, c,
	};
	return Matrix<T, 2, 2>(data);
}

template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const Vector<T, 3> &v) {
	Vector<T, 3> u = normalize(v);
	T s = sin(angle);
	T c = cos(angle);
	// note column-major order
	T data[] = {
		u[0] * u[0] + (1 - u[0] * u[0]) * c, u[0] * u[1] * (1 - c) + u[2] * s, u[0] * u[2] * (1 - c) - u[1] * s,
		u[0] * u[1] * (1 - c) - u[2] * s, u[1] * u[1] + (1 - u[1] * u[1]) * c, u[1] * u[2] * (1 - c) + u[0] * s,
		u[0] * u[2] * (1 - c) + u[1] * s, u[1] * u[2] * (1 - c) - u[0] * s, u[2] * u[2] + (1 - u[2] * u[2]) * c,
	};
	return Matrix<T, 3, 3>(data);
}

template <class T>
inline Matrix<T, 3, 3> rotationMatrix(const T &angle, const T &x, const T &y, const T &z){
	return rotationMatrix(angle, Vector<T, 3>(x, y, z));
}

template <class T, unsigned int M>
inline Matrix<T, M, M> scalingMatrix(const Vector<T, M> &v) {
	Matrix<T, M, M> ret(1);
	for (unsigned int i = 0; i < M; i++) {
		ret[i][i] = v[i];
	}
	return ret;
}

template <class T>
inline Matrix<T, 2, 2> scalingMatrix(const T &x, const T &y) {
	return scalingMatrix(Vector<T, 2>(x, y));
}

template <class T>
inline Matrix<T, 3, 3> scalingMatrix(const T &x, const T &y, const T &z) {
	return scalingMatrix(Vector<T, 3>(x, y, z));
}

template <class T>
inline Matrix<T, 4, 4> scalingMatrix(const T &x, const T &y, const T &z, const T &w) {
	return scalingMatrix(Vector<T, 4>(x, y, z, w));
}

template <class T, unsigned int M>
inline Matrix<T, M + 1, M + 1> translationMatrix(const Vector<T, M> &v) {
	Matrix<T, M + 1, M + 1> ret(1);
	for (unsigned int i = 0; i < M; i++) {
		ret[M][i] = v[i];
	}
	return ret;
}

template <class T>
inline Matrix<T, 3, 3> translationMatrix(const T &x, const T &y) {
	return translationMatrix(Vector<T, 2>(x, y));
}

template <class T>
inline Matrix<T, 4, 4> translationMatrix(const T &x, const T &y, const T &z) {
	return translationMatrix(Vector<T, 3>(x, y, z));
}

template <class T>
inline Matrix<T, 4, 4> perspectiveMatrix(const T &fovy, const T &aspect, const T &nearz, const T &farz) {
	T f = T(1) / tan(fovy / T(2));
	// note column-major order
	T data[] = {
		f / aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (farz + nearz) / (nearz - farz), T(-1),
		0, 0, T(2) * farz * nearz / (nearz - farz), 0
	};
	return Matrix<T, 4, 4>(data);
}

template <class T>
inline Matrix<T, 4, 4> frustumMatrix(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f) {
	T data[] = {
		T(2) * n / (r - l), 0, 0, 0,
		0, T(2) * n / (t - b), 0, 0,
		(r + l) / (r - l), (t + b) / (t - b), -((f + n) / (f - n)), T(-1),
		0, 0, T(-2) * f * n / (f - n), 0,
	};
	return Matrix<float, 4, 4>(data);
}

template <class T>
inline Matrix<T, 4, 4> orthoMatrix(const T &l, const T &r, const T &b, const T &t, const T &n, const T &f) {
	T data[] = {
		T(2) / (r - l), 0, 0, 0,
		0, T(2) / (t - b), 0, 0,
		0, 0, T(-2) / (f - n), 0,
		-((r + l) / (r - l)), -((t + b) / (t - b)), -((f + n) / (f - n)), T(1),
	};
	return Matrix<float, 4, 4>(data);
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator +=(const T &x) {
	Vector<T, M> col(x);
	for (unsigned int j = 0; j < N; j++) {
		(*this)[j] += col;
	}
	return *this;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> &Matrix<T, M, N>::operator -=(const T &x) {
	Vector<T, M> col(x);
	for (unsigned int j = 0; j < N; j++) {
		(*this)[j] -= col;
	}
	return *this;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const Matrix<T, M, N> &a, const T &b) {
	Matrix<T, M, N> ret(a);
	ret += b;
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator +(const T &a, const Matrix<T, M, N> &b) {
	Vector<T, M> col(a);
	Matrix<T, M, N> ret;
	for (unsigned int j = 0; j < N; j++) {
		ret[j] = col + b[j];
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const Matrix<T, M, N> &a, const T &b) {
	Matrix<T, M, N> ret(a);
	ret -= b;
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator -(const T &a, const Matrix<T, M, N> &b) {
	Vector<T, M> col(a);
	Matrix<T, M, N> ret;
	for (unsigned int j = 0; j < N; j++) {
		ret[j] = col - b[j];
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator /(const T &x, const Matrix<T, M, N> &a) {
	Vector<T, M> col(x);
	Matrix<T, M, N> ret;
	for (unsigned int j = 0; j < N; j++) {
		ret[j] = col / a[j];
	}
	return ret;
}

template <class T, unsigned int M, unsigned int N>
inline Matrix<T, M, N> operator *(const T &x, const Matrix<T, M, N> &a) {
	Vector<T, M> col(x);
	Matrix<T, M, N> ret;
	for (unsigned int j = 0; j < N; j++) {
		ret[j] = col / a[j];
	}
	return ret;
}



}

#endif //GLAM_MATRIX_H

