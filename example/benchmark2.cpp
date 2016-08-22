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

#include <iostream>
#include <random>
#include <vector>
#include <chrono>
//#include <glam/x86/vector_sse.h>
#include <glam/vector.h>

using namespace glam;

class TestSet {
private:
	size_t _size;
	std::vector<Vector<float, 2> > _pool;
	std::vector<unsigned int> _result;
	
public:
	template <class Generator>
	TestSet(size_t pool, Generator &gen) : _size(pool), _pool(pool), _result(pool) {
		for (auto it = _pool.begin(); it != _pool.end(); it++) {
			for (size_t i = 0; i < 2; i++) {
				(*it)[i] = std::generate_canonical<float, 24>(gen);
			}
		}
	}
	
	void operator() () {
		for (size_t i = 0; i < _size; i++) {
			_result[i] = packHalf2x16(_pool[i]);
		}
	}
};

int main(int argc, char **argv) {
	const size_t size = 100000000;
	std::random_device seed;
	std::default_random_engine rnd(seed());
	std::cout << "Generating random data for " << size << " vector pairs..." << std::endl;
	TestSet test(size, rnd);
	std::cout << "Running test..." << std::endl;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	test();
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "Duration: " << duration.count() << " seconds" << std::endl;
	std::cout << "Speed: " << (size / duration.count()) << " packs/second" << std::endl;
	return 0;
}

