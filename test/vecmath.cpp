/* Copyright (c) 2012, Gregor Riepl <onitake@gmail.com>
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

#include <cassert>
#include <iostream>
#include <glam/vector.h>
#include <glam/matrix.h>

int main(int argc, char **argv) {
	glam::ivec4 v(1, 12.345, 0.5, 30);
	std::cout << v << std::endl;
	glam::vec4 w(1, 12.345, 0.5, 30);
	std::cout << (w * glam::vec4(v)) << std::endl;
	glam::mat4 r = glam::orthoMatrix(-1.0f, 1.0f, -1.0f, 1.0f, 0.1f, 100.0f);
	std::cout << r << std::endl;
	std::cout << (r * w) << std::endl;
	glam::vec2 v1(3.0f, 10.0f);
	glam::vec2 v2(1.0f, 1.0f);
	float v3 = glam::dot(v1, v2);
	glam::vec2 v4 = v1 * glam::vec2(0.5f);
	glam::vec2 v5 = v1 + glam::vec2(1.0f, 2.0f);
	float x1 = v4[0];
	glam::vec4 v6;
	double values[] = { 10, 100.1, 42, 9.856, 19, 37, 37.3, 90, 101.85, 0 };
	glam::Vector<double, 10> v7(values);
	glam::Vector<double, 10> v8 = v7 * glam::Vector<double, 10>(73.123);
	std::cout << length(normalize(v8)) << std::endl;
	std::cout << v8 << std::endl;
	glam::Vector<double, 10> v9;
	double s1 = glam::dot(v8, v9);
	std::cout << s1 << std::endl;
	return 0;
}

