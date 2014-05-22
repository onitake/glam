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

#include <iostream>
#include <complex>
#include <vector>
#include <glam/matrix.h>

using namespace std;
using namespace glam;

const unsigned int N = 256;
const unsigned int S = N + N / 4;
typedef Vector<complex<double>, N> vec256cd;
typedef Matrix<complex<double>, N> mat256cd;
typedef Vector<double, N> vec256d;
typedef Vector<int, N> vec256i;

int main(int argc, char **argv) {
	// Calculate IQ separation coefficients
	/*vec256cd IQ;
	for (unsigned int n = 0; n < N; n++) {
		const int k = 1;
		const int p = 0;
		//IQ[n] = complex<double>(std::cos(k * 2 * M_PI * (n + p) / N), std::sin(k * 2 * M_PI * (n + p) / N));
		IQ[n] = complex<double>(1.0, n == 0 ? 0.0 : 1.0 / (2 * M_PI * n));
	}*/
	// Calculate discrete Hilbert transform operator
	mat256cd DHT;
	for (unsigned int k = 0; k < N; k++) {
		for (unsigned int n = 0; n < N; n++) {
			double d;
			if (k == n) {
				d = 1.0;
			} else {
				d = 0.0;
			}
			if (k & 1) {
				// k odd
				if (n & 1) {
					// n odd
					DHT[k][n] = complex<double>(d, 0);
				} else {
					// n even
					DHT[k][n] = complex<double>(d, 2 / ((k - n) * M_PI));
				}
			} else {
				// k even
				if (n & 1) {
					// n odd
					DHT[k][n] = complex<double>(d, 2 / ((k - n) * M_PI));
				} else {
					// n even
					DHT[k][n] = complex<double>(d, 0);
				}
			}
		}
	}
	//cout << "DHT=" << DHT << endl;
	// Calculate discrete Fourier transform operator
	complex<double> omega = exp(-2 * M_PI / N * complex<double>(0, 1));
	mat256cd DFT;
	for (unsigned int n = 0; n < N; n++) {
		for (unsigned int k = 0; k < N; k++) {
			DFT[n][k] = pow(omega, double(k * n));
		}
	}
	//cout << "DFT=" << DFT << endl;

	// Generate waveform
	//vector<double> V(N);
	vec256cd V;
	for (unsigned int n = 0; n < N; n++) {
		// Base frequency (as a multiple of 1 / (N * dt))
		const int k = 1;
		// Phase (as a multiple of 1 / N)
		const int p = 0;
		// Sine wave
		V[n] = cos(k * 2 * M_PI * (n + p) / N);
		// Square wave
		//V[n] = fract(double(k) * (n + p) / N) < 0.5 ? 0.0 : 1.0;
		// Positive sawtooth wave
		//V[n] = fract(double(k) * (n + p) / N);
		// Negative sawtooth wave
		//V[n] = 1.0 - fract(double(k) * (n + p) / N);
	}
	cout << "V=" << V << endl;

	// Determine analytic waveform
	//vec256d I(make_pair(&V[0], &V[N]));
	//vec256d Q(make_pair(&V[N / 4], &V[S]));
	//cout << "I=" << I << endl << "Q=" << Q << endl;
	//vec256cd Y = vec256cd(I) * vec256cd(complex<double>(1, 0)) + vec256cd(Q) * vec256cd(complex<double>(0, 1));
	vec256cd Y = DHT * V;
	cout << "Y=" << Y << endl;
	
	// Calculate the Fourier transform
	vec256cd X = DFT * Y;
	//cout << "X=" << X << endl;

	// Extract amplitude and phase angle
	vec256i L, A;
	for (unsigned int n = 0; n < N; n++) {
		L[n] = int(abs(X[n]));
		if (L[n] != 0) {
			A[n] = int(round(degrees(arg(X[n]))));
		}
	}
	// Print the result
	cout << "L=" << L << endl;
	cout << "A=" << A << endl;

	return 0;
}
