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
template <typename Type, size_t Height, size_t Width = Height>
class Matrix {
private:
	// The matrix is sorted in column-major order, as in OpenGL:
	// | a00 a01 a02 | => [ a00 a10 ] [ a01 a11 ] [ a02 a12 ] (index = Height * j + i)
	// | a10 a11 a12 |
	Vector<Type, Height> _a[Width];

	// Single element access by array index (inefficient implementation, use sparingly)
	Type &element(size_t index);
	inline const Type &element(size_t index) const;

	// Supporting initializers to avoid trouble from delegate constructors
	template <size_t Index>
	inline void initialize();
	template <size_t Index>
	inline void initialize(const Type &v);
	template <size_t Index, typename ForwardIterator>
	inline void initialize(std::pair<ForwardIterator, ForwardIterator> iterator);
	template <size_t Index>
	inline void initialize(const Type *v);
	template <size_t Index, typename TypeM, size_t HeightM, size_t WidthM>
	inline void initialize(const Matrix<TypeM, HeightM, WidthM> &m);
	template <size_t Index, typename... Args, size_t Size, typename Permutation>
	inline void initialize(const Vector<Type, Size, Permutation> &v, Args... args);
	template <size_t Index, typename... Args>
	inline void initialize(const Type &v, Args... args);
	
public:
	// Catch-all variadic constructor
	// The following constructor calls are supported:
	// Matrix(): Creates an empty matrix
	//  All elements are 0.
	// Matrix(const Type &d): Identity or scaling matrix
	//  The diagonal elements are all d, everything else is 0.
	// Matrix<ForwardIterator>(std::pair<ForwardIterator, ForwardIterator> iterator): Populate the matrix from an iterator pair
	//  Elements are filled up in order of occurrence, column by column, row by row. Transpose if you want the same layout as in a const array.
	//  std::pair is used to avoid call ambiguity problems. Call std::make_pair(begin, end) to constuct your iterator pair.
	// Matrix<TypeName, Height, Width>(const Matrix<TypeName, Height, Width> &m): Create a plain, cast or resized copy of another matrix
	//  Elements are filled up column by column, padded up with 0s. The diagonal is padded with 1s.
	// Matrix<AnyScalarOrVector...>(const AnyScalarOrVector &arg0, ...): Populate the matrix from a list of vectors and/or scalars
	//  Elements are filled up in order of occurrence, column by column, row by row.
	// Matrix(const Type *m): Populate the matrix from a constant array
	//  Elements are filled up in order of occurrence, column by column, row by row. Transpose if you want the same layout as in a const array.
	//  Note that this overload does not support range checking. Make sure that the input array has the correct size or use the ForwardIterator overload.
	template <typename... Args>
	inline Matrix(Args... args);
	
	// Get internal pointer, subject to memory alignment (padding between column vectors)
	inline const Type *internal() const;
	// Column vector access operator
	// Use double index operators (i.e. [j][i]) to access a single element
	inline Vector<Type, Height> &operator [](size_t j);
	// Const column vector access operator
	// Use double index operators (i.e. [j][i]) to access a single element
	inline const Vector<Type, Height> &operator [](size_t j) const;
	// Assignment operator
	inline Matrix<Type, Height, Width> &operator =(const Matrix<Type, Height, Width> &other);
	// Component-wise sum
	inline Matrix<Type, Height, Width> &operator +=(const Matrix<Type, Height, Width> &other);
	// Matrix-scalar sum
	inline Matrix<Type, Height, Width> &operator +=(const Type &x);
	// Component-wise difference
	inline 	Matrix<Type, Height, Width> &operator -=(const Matrix<Type, Height, Width> &other);
	// Matrix-scalar difference
	inline Matrix<Type, Height, Width> &operator -=(const Type &x);
	// Matrix division
	inline Matrix<Type, Height, Width> &operator /=(const Matrix<Type, Height, Width> &other);
	// Matrix-scalar division
	inline Matrix<Type, Height, Width> &operator /=(const Type &x);
	// Matrix-scalar product
	inline Matrix<Type, Height, Width> &operator *=(const Type &x);
};

// LU(P) decomposition of Height
// lower is a lower triangular matrix (L)
// upper is an upper triangular matrix (U)
// permutation is a permutation matrix (P)
// swaps is the number of swaps performed to obtain permutation, with respect to identity (S)
// The following equation must hold: PM = LU
template <typename Type, size_t Height, size_t Width>
struct LuDecomposition {
	Matrix<Type, Height, Width> lower;
	Matrix<Type, Height, Width> upper;
	Matrix<Type, Height, Width> permutation;
	size_t swaps;
	inline LuDecomposition(const Matrix<Type, Height, Width> &m);
};

// Component-wise matrix sum
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator +(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b);
// Component-wise matrix-scalar sum
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator +(const Matrix<Type, Height, Width> &a, const Type &b);
// Component-wise scalar-matrix sum
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator +(const Type &a, const Matrix<Type, Height, Width> &b);

// Component-wise matrix difference
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator -(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b);
// Component-wise matrix-scalar difference
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator -(const Matrix<Type, Height, Width> &a, const Type &b);
// Component-wise scalar-matrix difference
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator -(const Type &a, const Matrix<Type, Height, Width> &b);

// Component-wise matrix division
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator /(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b);
// Division of every matrix component by a scalar
template <typename Type, size_t Height, size_t Width, class U>
inline Matrix<Type, Height, Width> operator /(const Matrix<Type, Height, Width> &a, const U &x);
// Division of a scalar with every matrix component
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator /(const Type &x, const Matrix<Type, Height, Width> &a);

// Multiplication of every matrix component with a scalar
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator *(const Matrix<Type, Height, Width> &a, const Type &x);
// Multiplication of a scalar with every matrix component
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator *(const Type &x, const Matrix<Type, Height, Width> &a);

// Matrix-vector product
template <typename Type, size_t Height, size_t Width>
inline Vector<Type, Height> operator *(const Matrix<Type, Height, Width> &a, const Vector<Type, Width> &b);
// Vector-matrix product
template <typename Type, size_t Height, size_t Width>
inline Vector<Type, Height> operator *(const Vector<Type, Width> &a, const Matrix<Type, Height, Width> &b);

// Matrix product
template <typename Type, size_t Height, size_t Width, size_t WidthB>
inline Matrix<Type, Height, WidthB> operator *(const Matrix<Type, Height, Width> &a, const Matrix<Type, Width, WidthB> &b);

// Matrix comparison operator (exact comparison)
template <typename Type, size_t Height, size_t Width>
inline bool operator ==(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b);

// Component-wise matrix product
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> matrixCompMult(const Matrix<Type, Height, Width> &x, const Matrix<Type, Height, Width> &y);

// Cartesian vector product
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> outerProduct(const Vector<Type, Width> &c, const Vector<Type, Height> &r);

// Matrix transpose (columns and rows are swapped)
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Width, Height> transpose(const Matrix<Type, Height, Width> &m);

// Matrix determinant
template <typename Type, size_t Height, size_t Width>
inline Type determinant(const Matrix<Type, Height, Width> &m);

// Matrix inverse (through LU decomposition)
template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> inverse(const Matrix<Type, Height, Width> &m);

// Output stream operator
template <typename Type, size_t Height, size_t Width>
inline std::ostream &operator <<(std::ostream &o, const Matrix<Type, Height, Width> &m);

// Generate a rotation matrix
// Angle must be specified in radians.
// 2D linear transform
template <typename Type>
inline Matrix<Type, 2, 2> rotationMatrix(const Type &angle);
// 3D linear transform
template <typename Type>
inline Matrix<Type, 3, 3> rotationMatrix(const Type &angle, const Vector<Type, 3> &v);
template <typename Type>
inline Matrix<Type, 3, 3> rotationMatrix(const Type &angle, const Type &x, const Type &y, const Type &z);

// Generate a scaling matrix
// Linear transform
template <typename Type, size_t Height>
inline Matrix<Type, Height, Height> scalingMatrix(const Vector<Type, Height> &v);
// 2D linear transform
template <typename Type>
inline Matrix<Type, 2, 2> scalingMatrix(const Type &x, const Type &y);
// 3D linear transform
template <typename Type>
inline Matrix<Type, 3, 3> scalingMatrix(const Type &x, const Type &y, const Type &z);

// Generate a translation matrix
// Affine transform
template <typename Type, size_t Height>
inline Matrix<Type, Height + 1, Height + 1> translationMatrix(const Vector<Type, Height> &v);
// 2D affine transform
template <typename Type>
inline Matrix<Type, 3, 3> translationMatrix(const Type &x, const Type &y);
// 3D affine transform
template <typename Type>
inline Matrix<Type, 4, 4> translationMatrix(const Type &x, const Type &y, const Type &z);

// Generate a perspective projection matrix
// Specification through field-of-view angle, aspect ratio and near/far planes.
// Angle must be specified in radians.
template <typename Type>
inline Matrix<Type, 4, 4> perspectiveMatrix(const Type &fovy, const Type &aspect, const Type &nearz, const Type &farz);
// Specification through view frustum edges (left, right, bottom, top, near, far)
template <typename Type>
inline Matrix<Type, 4, 4> frustumMatrix(const Type &l, const Type &r, const Type &b, const Type &t, const Type &n, const Type &f);

// Generate an orthographic projection matrix
// Specification through view cube edges (left, right, bottom, top, near, far).
template <typename Type>
inline Matrix<Type, 4, 4> orthoMatrix(const Type &l, const Type &r, const Type &b, const Type &t, const Type &n, const Type &f);

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

template <typename Type, size_t Height, size_t Width>
inline Type &Matrix<Type, Height, Width>::element(size_t index) {
	return (*this)[index / Height][index % Height];
}

template <typename Type, size_t Height, size_t Width>
inline const Type &Matrix<Type, Height, Width>::element(size_t index) const {
	return (*this)[index / Height][index % Height];
}

template <typename Type, size_t Height, size_t Width>
template <size_t Index>
inline void Matrix<Type, Height, Width>::initialize() { }

template <typename Type, size_t Height, size_t Width>
template <size_t Index, typename... Args>
inline void Matrix<Type, Height, Width>::initialize(const Type &v, Args... args) {
	// The condition is < and not <= because there is a dedicated terminal initializer for a single scalar
	static_assert(Index + 1 < Height * Width, "Too many initializers");
	element(Index) = v;
	initialize<Index + 1>(args...);
}

template <typename Type, size_t Height, size_t Width>
template <size_t Index, typename... Args, size_t Size, typename Permutation>
inline void Matrix<Type, Height, Width>::initialize(const Vector<Type, Size, Permutation> &v, Args... args) {
	static_assert(Index + Permutation::Elements <= Height * Width, "Too many initializers");
	for (size_t p = 0; p < Permutation::Elements; p++) {
		element(Index + p) = v[p];
	}
	initialize<Index + Permutation::Elements>(args...);
}

template <typename Type, size_t Height, size_t Width>
template <size_t Index>
inline void Matrix<Type, Height, Width>::initialize(const Type &v) {
	static_assert(Index + 1 <= Height * Width, "Too many initializers");
	if (Index == 0) {
		for (size_t i = 0; i < Height; i++) {
			for (size_t j = 0; j < Width; j++) {
				if (i == j) {
					// Diagonal element, fill with value
					(*this)[j][i] = v;
				} else {
					// Anything else, fill with 0
					(*this)[j][i] = Type(0);
				}
			}
		}
	} else {
		element(Index) = v;
	}
}

template <typename Type, size_t Height, size_t Width>
template <size_t Index, typename ForwardIterator>
inline void Matrix<Type, Height, Width>::initialize(std::pair<ForwardIterator, ForwardIterator> iterator) {
	static_assert(Index == 0, "No other arguments are allowed when initializing a matrix from an iterator");
	size_t index = Index;
	for (ForwardIterator it = iterator.first; it != iterator.second; it++, index++) {
		element(index) = *it;
	}
}

template <typename Type, size_t Height, size_t Width>
template <size_t Index>
inline void Matrix<Type, Height, Width>::initialize(const Type *v) {
	static_assert(Index == 0, "No other arguments are allowed when initializing a matrix from an array");
	for (size_t p = 0; p < Height * Width; p++) {
		element(Index + p) = v[p];
	}
}

template <typename Type, size_t Height, size_t Width>
template <size_t Index, typename TypeM, size_t HeightM, size_t WidthM>
inline void Matrix<Type, Height, Width>::initialize(const Matrix<TypeM, HeightM, WidthM> &m) {
	static_assert(Index == 0, "No other arguments are allowed when copy-constructing a matrix");
	for (size_t i = 0; i < Height; i++) {
		for (size_t j = 0; j < Width; j++) {
			if (i < HeightM && j < WidthM) {
				// Component available, copy
				(*this)[j][i] = Type(m[j][i]);
			} else {
				// Component not available, fill
				if (i == j) {
					// Diagonal element, fill with identity
					(*this)[j][i] = Type(1);
				} else {
					// Anything else, fill with 0
					(*this)[j][i] = Type(0);
				}
			}
		}
	}
}

template <typename Type, size_t Height, size_t Width>
template <typename... Args>
inline Matrix<Type, Height, Width>::Matrix(Args... args) {
	initialize<0>(args...);
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator =(const Matrix<Type, Height, Width> &other) {
	for (size_t j = 0; j < Width; j++) {
		(*this)[j] = other[j];
	}
	return *this;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator +=(const Matrix<Type, Height, Width> &other) {
	for (size_t j = 0; j < Width; j++) {
		(*this)[j] += other[j];
	}
	return *this;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator +(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b) {
	Matrix<Type, Height, Width> ret = a;
	ret += b;
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator -=(const Matrix<Type, Height, Width> &other) {
	for (size_t j = 0; j < Width; j++) {
		(*this)[j] -= other[j];
	}
	return *this;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator -(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b) {
	Matrix<Type, Height, Width> ret = a;
	ret -= b;
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator *=(const Type &x) {
	for (size_t i = 0; i < Height * Width; i++) {
		(*this)[i] *= x;
	}
	return *this;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator *(const Matrix<Type, Height, Width> &a, const Type &x) {
	Matrix<Type, Height, Width> ret = a;
	ret *= x;
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator /=(const Type &x) {
	*this *= 1 / x;
	return *this;
}

template <typename Type, size_t Height, size_t Width, class U>
inline Matrix<Type, Height, Width> operator /(const Matrix<Type, Height, Width> &a, const U &x) {
	Matrix<Type, Height, Width> ret = a;
	ret /= x;
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline std::ostream &operator <<(std::ostream &o, const Matrix<Type, Height, Width> &m) {
	o << "{";
	for (size_t i = 0; i < Height; i++) {
		o << "\n\t";
		for (size_t j = 0; j < Width; j++) {
			o << m[j][i] << "\t";
		}
	}
	o << "\n}";
	return o;
}

template <typename Type, size_t Height, size_t Width>
inline bool operator ==(const Matrix<Type, Height, Width> &a, const Matrix<Type, Height, Width> &b) {
	for (size_t j = 0; j < Width; j++) {
		if (a[j] != b[j]) {
			return false;
		}
	}
	/*for (size_t i = 0; i < Height; i++) {
		for (size_t j = 0; j < Width; j++) {
			if (a[j][i] != b[j][i]) {
				return false;
			}
		}
	}*/
	return true;
}

template <typename Type, size_t Height, size_t Width, size_t WidthB>
inline Matrix<Type, Height, WidthB> operator *(const Matrix<Type, Height, Width> &a, const Matrix<Type, Width, WidthB> &b) {
	//      4 x 3               3 x 4              4 x 4
	// | a00 a01 a02 |   | b00 b01 b02 b03 |   | a00*b00+a01*b10+a02*b20 a00*b01+a01*b11+a02*b21 ... |
	// | a10 a11 a12 | * | b10 b11 b12 b33 | = | a10*b00+a11*b10+a12*b20 ...                         |
	// | a20 a21 a22 |   | b20 b21 b22 b33 |   | ...                                                 |
	// | a30 a31 a32 |                         | ...                                                 |
	Matrix<Type, Height, WidthB> ret;
	Matrix<Type, Width, Height> t = transpose(a);
	for (size_t i = 0; i < Height; i++) {
		for (size_t j = 0; j < WidthB; j++) {
			ret[j][i] = dot(t[i], b[j]);
		}
	}
	/*for (size_t i = 0; i < Height; i++) {
		for (size_t j = 0; j < P; j++) {
			ret[j][i] = 0;
			for (size_t k = 0; k < Width; k++) {
				ret[j][i] += a[k][i] * b[j][k];
			}
		}
	}*/
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Vector<Type, Height> operator *(const Matrix<Type, Height, Width> &a, const Vector<Type, Width> &b) {
	// Transpose the matrix
	Vector<Type, Height> ret;
	Matrix<Type, Width, Height> t = transpose(a);
	for (size_t i = 0; i < Height; i++) {
		ret[i] = dot(t[i], b);
	}
	/*for (size_t i = 0; i < Height; i++) {
		ret[i] = 0;
		for (size_t k = 0; k < Width; k++) {
			ret[i] += a[k][i] * b[k];
		}
	}*/
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline const Type *Matrix<Type, Height, Width>::internal() const {
	return &(*this)[0][0];
}

template <typename Type, size_t Height, size_t Width>
inline Vector<Type, Height> &Matrix<Type, Height, Width>::operator [](size_t j) {
#ifdef GLAM_RANGE_CHECKS
	if (j < Width) {
		return _a[j];
	}
	throw DimensionOutOfRangeException("Matrix column index out of range", Width, j);
#else
	return _a[j];
#endif
}

template <typename Type, size_t Height, size_t Width>
inline const Vector<Type, Height> &Matrix<Type, Height, Width>::operator [](size_t j) const {
#ifdef GLAM_RANGE_CHECKS
	if (j < Width) {
		return _a[j];
	}
	throw DimensionOutOfRangeException("Matrix column index out of range", Width, j);
#else
	return _a[j];
#endif
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> matrixCompMult(const Matrix<Type, Height, Width> &x, const Matrix<Type, Height, Width> &y) {
	Matrix<Type, Height, Width> ret;
	for (size_t i = 0; i < Height; i++) {
		for (size_t j = 0; j < Width; j++) {
			ret[j][i] = x[j][i] * y[j][i];
		}
	}
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> outerProduct(const Vector<Type, Height> &c, const Vector<Type, Width> &r) {
	return Matrix<Type, Height, 1>(c) * Matrix<Type, 1, Width>(r);
}

template <typename Type, size_t Height, size_t Width>
inline LuDecomposition<Type, Height, Width>::LuDecomposition(const Matrix<Type, Height, Width> &m) : lower(1), upper(m), permutation(1), swaps(0) {
	static_assert(Height == Width, "Matrix is not square");
	// Initializer list: prepare result (lower and upper triangular matrices, permutation, reset the swap counter)
	// Loop through all rows except the last (which will simply become 0 0 .. 0 1 in L)
	for (size_t n = 0; n < Height - 1; n++) {
		// Prepare the intermediate lower matrix (one column filled)
		Matrix<Type, Height, Width> ln(1);
		// Check if a swap is needed
		size_t s = Height;
		for (size_t i = n; i < Height && s == Height; i++) {
			if (upper[i][n] != 0) {
				s = i;
			}
		}
		if (s == Height) {
			// All coefficients of column n are 0
			throw NonInvertibleMatrixException("Singular matrix, can't find solution for inverse");
		} else if (s != n) {
			// (n, n) is 0, swap row n with row s (which has a non-zero coefficient)
			// Also swap the corresponding rows in the permutation matrix
			for (size_t j = 0; j < Width; j++) {
				std::swap(upper[j][n], upper[j][s]);
				std::swap(permutation[j][n], permutation[j][s]);
			}
			// Increment the swap counter
			swaps++;
		}
		// Calculate lower matrix coefficients for this column
		for (size_t i = n + 1; i < Height; i++) {
			Type v = upper[n][i] / upper[n][n];
			// And assign them to the appropriate Ln
			ln[n][i] = -v;
			// and L fields
			lower[n][i] = v;
		}
		// Apply Ln to U, this cancels out all coefficients under U(n,n)
		upper = ln * upper;
	}
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Width, Height> transpose(const Matrix<Type, Height, Width> &m) {
	Matrix<Type, Width, Height> ret;
	for (size_t i = 0; i < Height; i++) {
		for (size_t j = 0; j < Width; j++) {
			ret[i][j] = m[j][i];
		}
	}
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Type determinant(const Matrix<Type, Height, Width> &m) {
	static_assert(Height == Width && Height > 0, "Matrix is not square");
	// Decompose into triangular parts
	try {
		LuDecomposition<Type, Height, Width> lups(m);
		// det(A) = det(P^-1) * det(L) * det(U) = (-1)^S * (l11 * l22 * ...) * (u11 * u22 * ...)
		// Since the diagonal of L is all 1s, its determinant is also 1, so
		// det(A) = (-1)^S * (u11 * u22 * ...)
		// Calculate d according to (-1)^S here, S = number of exchanges in P^-1
		Type d;
		if (lups.swaps & 1) {
			d = Type(-1);
		} else {
			d = Type(1);
		}
		for (size_t i = 0; i < Height; i++) {
			d *= lups.upper[i][i];
		}
		return d;
	} catch (NonInvertibleMatrixException &ex) {
		return Type(0);
	}
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> inverse(const Matrix<Type, Height, Width> &m) {
	static_assert(Height == Width && Height > 0, "Matrix is not square");
	// Decompose into triangular parts
	LuDecomposition<Type, Height, Width> lups(m);
	// Start with AX = I, where Index is the inverse of A
	// A = P^-1LU, so P^-1LUX = I, and thus LUX = P (with P = I if no row swapping was needed)
	// Substitute UX = Y, yielding LY = P
	// Calculate Y through forward substitution of L
	Matrix<Type, Height, Width> y;
	// Process each column independently
	for (size_t j = 0; j < Width; j++) {
		// Substitute each row from the previous rows
		for (size_t i = 0; i < Height; i++) {
			// Fetch the starting value from P
			Type diff = lups.permutation[j][i];
			// Subtract the product of the values from the previous rows (from the top) times the corresponding values in L
			for (size_t k = 0; k < i; k++) {
				diff -= lups.lower[k][i] * y[j][k];
			}
			// No division necessary, the diagonal is always 1
			y[j][i] = diff;
		}
	}
	// Calculate Index from UX = Y through backward substitution of U
	Matrix<Type, Height, Width> x;
	// Process each column independently
	for (size_t j = 0; j < Width; j++) {
		// Substitute each row from the previous rows, work backwards
		for (size_t i = Height; i-- > 0;) {
			// Fetch the starting value from Y
			Type diff = y[j][i];
			// Subtract the product of the values from the previous rows (from the bottom) times the corresponding values in U
			for (size_t k = i + 1; k < Height; k++) {
				diff -= lups.upper[k][i] * x[j][k];
			}
			// Divide by the diagonal coefficient
			x[j][i] = diff / lups.upper[i][i];
		}
	}
	// Index is the inverse of A
	return x;
}

template <typename Type>
inline Matrix<Type, 2, 2> rotationMatrix(const Type &angle) {
	Type s = sin(angle);
	Type c = cos(angle);
	// note column-major order
	Type data[] = {
		c, s,
		-s, c,
	};
	return Matrix<Type, 2, 2>(data);
}

template <typename Type>
inline Matrix<Type, 3, 3> rotationMatrix(const Type &angle, const Vector<Type, 3> &v) {
	Vector<Type, 3> u = normalize(v);
	Type s = sin(angle);
	Type c = cos(angle);
	// note column-major order
	Type data[] = {
		u[0] * u[0] + (1 - u[0] * u[0]) * c, u[0] * u[1] * (1 - c) + u[2] * s, u[0] * u[2] * (1 - c) - u[1] * s,
		u[0] * u[1] * (1 - c) - u[2] * s, u[1] * u[1] + (1 - u[1] * u[1]) * c, u[1] * u[2] * (1 - c) + u[0] * s,
		u[0] * u[2] * (1 - c) + u[1] * s, u[1] * u[2] * (1 - c) - u[0] * s, u[2] * u[2] + (1 - u[2] * u[2]) * c,
	};
	return Matrix<Type, 3, 3>(data);
}

template <typename Type>
inline Matrix<Type, 3, 3> rotationMatrix(const Type &angle, const Type &x, const Type &y, const Type &z){
	return rotationMatrix(angle, Vector<Type, 3>(x, y, z));
}

template <typename Type, size_t Height>
inline Matrix<Type, Height, Height> scalingMatrix(const Vector<Type, Height> &v) {
	Matrix<Type, Height, Height> ret(1);
	for (size_t i = 0; i < Height; i++) {
		ret[i][i] = v[i];
	}
	return ret;
}

template <typename Type>
inline Matrix<Type, 2, 2> scalingMatrix(const Type &x, const Type &y) {
	return scalingMatrix(Vector<Type, 2>(x, y));
}

template <typename Type>
inline Matrix<Type, 3, 3> scalingMatrix(const Type &x, const Type &y, const Type &z) {
	return scalingMatrix(Vector<Type, 3>(x, y, z));
}

template <typename Type>
inline Matrix<Type, 4, 4> scalingMatrix(const Type &x, const Type &y, const Type &z, const Type &w) {
	return scalingMatrix(Vector<Type, 4>(x, y, z, w));
}

template <typename Type, size_t Height>
inline Matrix<Type, Height + 1, Height + 1> translationMatrix(const Vector<Type, Height> &v) {
	Matrix<Type, Height + 1, Height + 1> ret(1);
	for (size_t i = 0; i < Height; i++) {
		ret[Height][i] = v[i];
	}
	return ret;
}

template <typename Type>
inline Matrix<Type, 3, 3> translationMatrix(const Type &x, const Type &y) {
	return translationMatrix(Vector<Type, 2>(x, y));
}

template <typename Type>
inline Matrix<Type, 4, 4> translationMatrix(const Type &x, const Type &y, const Type &z) {
	return translationMatrix(Vector<Type, 3>(x, y, z));
}

template <typename Type>
inline Matrix<Type, 4, 4> perspectiveMatrix(const Type &fovy, const Type &aspect, const Type &nearz, const Type &farz) {
	Type f = Type(1) / tan(fovy / Type(2));
	// note column-major order
	Type data[] = {
		f / aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (farz + nearz) / (nearz - farz), Type(-1),
		0, 0, Type(2) * farz * nearz / (nearz - farz), 0
	};
	return Matrix<Type, 4, 4>(data);
}

template <typename Type>
inline Matrix<Type, 4, 4> frustumMatrix(const Type &l, const Type &r, const Type &b, const Type &t, const Type &n, const Type &f) {
	Type data[] = {
		Type(2) * n / (r - l), 0, 0, 0,
		0, Type(2) * n / (t - b), 0, 0,
		(r + l) / (r - l), (t + b) / (t - b), -((f + n) / (f - n)), Type(-1),
		0, 0, Type(-2) * f * n / (f - n), 0,
	};
	return Matrix<float, 4, 4>(data);
}

template <typename Type>
inline Matrix<Type, 4, 4> orthoMatrix(const Type &l, const Type &r, const Type &b, const Type &t, const Type &n, const Type &f) {
	Type data[] = {
		Type(2) / (r - l), 0, 0, 0,
		0, Type(2) / (t - b), 0, 0,
		0, 0, Type(-2) / (f - n), 0,
		-((r + l) / (r - l)), -((t + b) / (t - b)), -((f + n) / (f - n)), Type(1),
	};
	return Matrix<float, 4, 4>(data);
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator +=(const Type &x) {
	Vector<Type, Height> col(x);
	for (size_t j = 0; j < Width; j++) {
		(*this)[j] += col;
	}
	return *this;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> &Matrix<Type, Height, Width>::operator -=(const Type &x) {
	Vector<Type, Height> col(x);
	for (size_t j = 0; j < Width; j++) {
		(*this)[j] -= col;
	}
	return *this;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator +(const Matrix<Type, Height, Width> &a, const Type &b) {
	Matrix<Type, Height, Width> ret(a);
	ret += b;
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator +(const Type &a, const Matrix<Type, Height, Width> &b) {
	Vector<Type, Height> col(a);
	Matrix<Type, Height, Width> ret;
	for (size_t j = 0; j < Width; j++) {
		ret[j] = col + b[j];
	}
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator -(const Matrix<Type, Height, Width> &a, const Type &b) {
	Matrix<Type, Height, Width> ret(a);
	ret -= b;
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator -(const Type &a, const Matrix<Type, Height, Width> &b) {
	Vector<Type, Height> col(a);
	Matrix<Type, Height, Width> ret;
	for (size_t j = 0; j < Width; j++) {
		ret[j] = col - b[j];
	}
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator /(const Type &x, const Matrix<Type, Height, Width> &a) {
	Vector<Type, Height> col(x);
	Matrix<Type, Height, Width> ret;
	for (size_t j = 0; j < Width; j++) {
		ret[j] = col / a[j];
	}
	return ret;
}

template <typename Type, size_t Height, size_t Width>
inline Matrix<Type, Height, Width> operator *(const Type &x, const Matrix<Type, Height, Width> &a) {
	Vector<Type, Height> col(x);
	Matrix<Type, Height, Width> ret;
	for (size_t j = 0; j < Width; j++) {
		ret[j] = col / a[j];
	}
	return ret;
}



}

#endif //GLAM_MATRIX_H

