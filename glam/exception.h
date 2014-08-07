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

#ifndef GLAM_EXCEPTION_H
#define GLAM_EXCEPTION_H

#include <exception>
#include <string>
#include <sstream>
#include <ostream>

namespace glam {

using std::exception;
using std::string;
using std::stringstream;
using std::ostream;

class Exception : public exception {
private:
	string _cause;

protected:
	virtual string &cause() throw() {
		return _cause;
	}

public:
	Exception() throw() { }
	Exception(const string &cause) throw() : _cause(cause) { }
	virtual ~Exception() throw() { }
	virtual const char* what() const throw() {
		return _cause.c_str();
	}
	virtual const string &cause() const throw() {
		return _cause;
	}
	virtual string type() const {
		return "Exception";
	}
};

inline ostream &operator <<(ostream &os, const Exception &e);

class DimensionOutOfRangeException : public Exception {
public:
	DimensionOutOfRangeException() : Exception() { }
	DimensionOutOfRangeException(const string &cause) : Exception(cause) { }
	DimensionOutOfRangeException(const string &cause, unsigned int max, unsigned int index) {
		stringstream cat;
		cat << cause << " (" << index << " > " << max << ")";
		this->cause() = cat.str();
	}
	string type() const {
		return "DimensionOutOfRangeException";
	}
};

class InvalidArgumentException : public Exception {
public:
	InvalidArgumentException() : Exception() { }
	InvalidArgumentException(const string &cause) : Exception(cause) { }
	string type() const {
		return "InvalidArgumentException";
	}
};

class NonInvertibleMatrixException : public Exception {
public:
	NonInvertibleMatrixException() : Exception() { }
	NonInvertibleMatrixException(const string &cause) : Exception(cause) { }
	string type() const {
		return "NonInvertibleMatrixException";
	}
};

inline ostream &operator <<(ostream &os, const Exception &e) {
	os << "[" << e.type() << "] " << e.cause();
	return os;
}

}

#endif //GLAM_EXCEPTION_H

