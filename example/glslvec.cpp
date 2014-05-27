/* Copyright (c) 2013, Gregor Riepl <onitake@gmail.com>
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
#include <glam/vector.h>
#include <glam/matrix.h>

using namespace std;
using namespace glam;

int main(int argc, char **argv) {
	vec2 a(1.0f, 0.5f);
	vec2 b = normalize(a);
	vec3 i;
	cout << "a=" << a << " b.xyxy=" << b.xyxy() << " i=" << i << endl;
	vec3 N(1, 2, 3);
	vec3 T(0.5, 1, 0);
	cout << "N·T=" << dot(N, T) << " NxT=" << cross(N, T) << " NxT·N=" << dot(cross(N, T), N) << endl;
	static const float _M[] = {
		-1, 0, 0, 0,
		0, -1, 0, 0,
		0, 0, -1, 0,
		0, 0, 0, 1,
	};
	mat4 I(1);
	mat4 M(_M);
	mat4 O(std::make_pair(_M, &_M[16]));
	cout << "M" << (M == O ? "==" : "!=") << "O" << endl;
	cout << "M·I=" << (M * I) << endl;
	cout << "M·I·vec4(N)=" << (M * I * vec4(N, 1)) << endl;
	Matrix<int, 2, 2> M3A(10, 20, 30, 40);
	Matrix<int, 2, 1> M3B(ivec2(10, 20));
	cout << "M3A·M3B=" << (M3A * M3B) << endl;
	return 0;
}
