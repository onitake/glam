/*
 * GLAM - GLSL Linear Algebra Math Library
 * 
 * Copyright (c) 2013, Gregor Riepl <onitake@gmail.com>
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
#include <fstream>
#include <vector>
#include <glam/vector.h>
#include <glam/matrix.h>

#include "sky_d.h"
#include "sky_n.h"
#include "teapot.h"

using namespace glam;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

struct sampler2D {
	ivec2 size;
	vec4 *samples;
};

struct Uniform {
	vec4 Kd;
	vec4 Ks;
	vec3 Ka;
	vec3 lightColor;
	vec3 ambientColor;
	sampler2D texture_d;
	sampler2D texture_n;
	mat4 projectionModelViewMatrix;
	mat3 normalMatrix;
	vec3 objectLightPosition;
};

struct Attribute {
	vec3 position;
	vec2 texc;
	vec3 normal;
	vec4 tangent;
};

struct Varying {
	vec4 position;
	vec2 texcOut;
	vec3 tangentLightOut;
	vec3 tangentPositionOut;
};

struct Output {
	vec4 fragColor;
};

struct Triangle {
	Varying a;
	Varying b;
	Varying c;
};

sampler2D header2sampler(unsigned int width, unsigned int height, const unsigned char *data) {
	sampler2D ret;
	ret.size = ivec2(int(width), int(height));
	unsigned int size = ret.size.x * ret.size.y;
	ret.samples = new vec4[size];
	for (unsigned int i = 0; i < size; i++) {
		ret.samples[i] = vec4(data[i * 3 + 0] / 255.0f, data[i * 3 + 1] / 255.0f, data[i * 3 + 2] / 255.0f, 1.0f);
	}
	return ret;
}

void bezier2mesh(const vector<unsigned int> &patches, const vector<float> &vertices, const uvec2 &subdivisions, vector<unsigned int> &indices, vector<Attribute> &attributes) {
	unsigned int offset = attributes.size();
	for (size_t v = 0; v < vertices.size() - 2; v += 3) {
		Attribute attr;
		attr.position = vec3(vertices[v], vertices[v + 1], vertices[v + 2]);
		attributes.push_back(attr);
	}
	vector<unsigned int> nscale(attributes.size());
	for (size_t b = 0; b < patches.size() - 15; b += 16) {
		for (size_t j = 0; j < 3; j++) {
			for (size_t i = 0; i < 3; i++) {
				uvec4 quad(b + j * 4 + i);
				quad.t += 1;
				quad.p += 4;
				quad.q += 5;
				uvec4 index(patches[quad.s], patches[quad.t], patches[quad.p], patches[quad.q]);
				/*if (index.s == index.t || index.s == index.p || index.s == index.q || index.t == index.p || index.t == index.q || index.p == index.q) {
					cout << "Coordinate is refered to twice: " << index << endl;
				}*/
				uvec4 oindex = index + offset;
				indices.push_back(oindex.s);
				indices.push_back(oindex.p);
				indices.push_back(oindex.t);
				indices.push_back(oindex.p);
				indices.push_back(oindex.q);
				indices.push_back(oindex.t);
				vec3 s = attributes[oindex.s].position;
				vec3 t = attributes[oindex.s].position;
				vec3 p = attributes[oindex.s].position;
				vec3 q = attributes[oindex.s].position;
				attributes[oindex.s].texc = vec2(i + 0.0f, j + 0.0f) / 4.0f;
				attributes[oindex.t].texc = vec2(i + 1.0f, j + 0.0f) / 4.0f;
				attributes[oindex.p].texc = vec2(i + 0.0f, j + 1.0f) / 4.0f;
				attributes[oindex.q].texc = vec2(i + 1.0f, j + 1.0f) / 4.0f;
				vec3 stpnormal = normalize(cross(t - s, p - s));
				vec3 stptangent = normalize(t - s);
				vec3 stpbitangent = normalize(p - s);
				float stpsign = sign(dot(cross(stpnormal, stptangent), stpbitangent));
				vec3 qptnormal = normalize(cross(q - p, t - p));
				vec3 qpttangent = normalize(q - p);
				vec3 qptbitangent = normalize(q - t);
				float qptsign = sign(dot(cross(qptnormal, qpttangent), qptbitangent));
				attributes[oindex.s].normal += stpnormal;
				attributes[oindex.s].tangent.xyz += stptangent;
				attributes[oindex.s].tangent.w *= stpsign;
				nscale[index.s]++;
				attributes[oindex.t].normal += stpnormal;
				attributes[oindex.t].tangent.xyz += stptangent;
				attributes[oindex.t].tangent.w *= stpsign;
				nscale[index.t]++;
				attributes[oindex.p].normal += stpnormal;
				attributes[oindex.p].tangent.xyz += stptangent;
				attributes[oindex.p].tangent.w *= stpsign;
				nscale[index.p]++;
				attributes[oindex.q].normal += qptnormal;
				attributes[oindex.q].tangent.xyz += qpttangent;
				attributes[oindex.q].tangent.w *= qptsign;
				nscale[index.q]++;
				attributes[oindex.p].normal += qptnormal;
				attributes[oindex.p].tangent.xyz += qpttangent;
				attributes[oindex.p].tangent.w *= qptsign;
				nscale[index.p]++;
				attributes[oindex.t].normal += qptnormal;
				attributes[oindex.t].tangent.xyz += qpttangent;
				attributes[oindex.t].tangent.w *= qptsign;
				nscale[index.t]++;
			}
		}
	}
	for (size_t n = 0; n < nscale.size(); n++) {
		Attribute &attr = attributes[n + offset];
		float scale = nscale[n] > 0 ? 1.0f / nscale[n] : 1.0f;
		attr.normal *= scale;
		attr.tangent *= scale;
		attr.tangent.w /= length(cross(attr.normal, attr.tangent.xyz));
	}
	/*Attribute vertices[] = {
		{ .position = vec3(-1.0f, 1.0f, 0.0f), .normal = normalize(vec3(-1.0f, 1.0f, 1.0f)), .texc = vec2(0.0f, 0.0f), .tangent = vec4() },
		{ .position = vec3(1.0f, 1.0f, 0.0f), .normal = normalize(vec3(1.0f, 1.0f, 1.0f)), .texc = vec2(1.0f, 0.0f), .tangent = vec4() },
		{ .position = vec3(1.0f, -1.0f, 0.0f), .normal = normalize(vec3(1.0f, -1.0f, 1.0f)), .texc = vec2(1.0f, 1.0f), .tangent = vec4() },
		{ .position = vec3(-1.0f, -1.0f, 0.0f), .normal = normalize(vec3(-1.0f, -1.0f, 1.0f)), .texc = vec2(0.0f, 1.0f), .tangent = vec4() },
	};
	unsigned int indices[] = { 0, 1, 2, 2, 3, 0 };*/
}

vec4 blend(const vec4 &tl, const vec4 &tr, const vec4 &bl, const vec4 &br, const vec2 &position) {
	vec4 t = mix(tl, tr, position.x);
	vec4 b = mix(bl, br, position.x);
	return mix(t, b, position.y);
}

Varying blend(const Varying &v0, const Varying &v1, const Varying &v2, const vec3 &bias) {
	Varying ret;
	ret.position = v0.position * bias.x + v1.position * bias.y + v2.position * bias.z;
	ret.texcOut = v0.texcOut * bias.x + v1.texcOut * bias.y + v2.texcOut * bias.z;
	ret.tangentLightOut = v0.tangentLightOut * bias.x + v1.tangentLightOut * bias.y + v2.tangentLightOut * bias.z;
	ret.tangentPositionOut = v0.tangentPositionOut * bias.x + v1.tangentPositionOut * bias.y + v2.tangentPositionOut * bias.z;
	return ret;
}

void texelSet(sampler2D &s, const ivec2 &p, int lod, const vec4 &c) {
	if (p.x < 0 || p.x >= s.size.x || p.y < 0 || p.y >= s.size.x) {
		cout << "Invalid sampler coordinates: " << p << " maximum: " << s.size << endl;
	} else {
		s.samples[p.y * s.size.x + p.x] = c;
	}
}

vec4 texelFetch(const sampler2D &s, const ivec2 &p, int lod) {
	if (p.x < 0 || p.x >= s.size.x || p.y < 0 || p.y >= s.size.x) {
		return vec4();
	}
	return s.samples[p.y * s.size.x + p.x];
}

vec4 texture(const sampler2D &s, const vec2 &p) {
	vec2 a = p * vec2(s.size);
	vec2 o;
	vec2 f = modf(a, o);
	vec4 t0 = texelFetch(s, ivec2(o) + ivec2(0, 0), 0);
	vec4 t1 = texelFetch(s, ivec2(o) + ivec2(1, 0), 0);
	vec4 t2 = texelFetch(s, ivec2(o) + ivec2(0, 1), 0);
	vec4 t3 = texelFetch(s, ivec2(o) + ivec2(1, 1), 0);
	return blend(t0, t1, t2, t3, f);
}

Varying vertexShader(const Uniform &u, const Attribute &a) {
	mat3 tangentMatrix = transpose(mat3(a.tangent.xyz, cross(a.normal, a.tangent.xyz) * a.tangent.w, a.normal));
	Varying v;
	v.tangentLightOut = tangentMatrix * (u.objectLightPosition - a.position);
	v.tangentPositionOut = tangentMatrix * a.position;
	v.position = u.projectionModelViewMatrix * vec4(a.position, 1.0f / a.position.z);
	v.texcOut = a.texc;
	v.position = vec4(a.position, 1.0f / a.position.z);
	return v;
}

Output fragmentShader(const Uniform &u, const Varying &v) {
	vec4 normal = texture(u.texture_n, v.texcOut);
	vec4 material = texture(u.texture_d, v.texcOut);
	vec3 N = normal.xyz * 2.0f - 1.0f;
	vec3 L = normalize(v.tangentLightOut);
	vec3 V = normalize(v.tangentPositionOut);
	vec3 R = (L + V) * 0.5f;
	vec3 ambient = u.ambientColor * u.Ka;
	vec3 diffuse = u.lightColor * u.Kd.xyz * clamp(dot(L, N), 0.0f, 1.0f);
	vec3 specular = u.lightColor * u.Ks.xyz * pow(clamp(dot(R, N), 0.0f, 1.0f), u.Ks.w);
	Output o;
	o.fragColor = vec4((ambient + diffuse) * material.xyz + specular, u.Kd.w * material.w);
	vec2 color = v.position.xy * vec2(0.5) + vec2(0.5);
	o.fragColor = vec4(color, 0.0f, 1.0f);
	return o;
}

vector<Varying> vertexStage(const Uniform &uniform, const vector<Attribute> &attributes) {
	vector<Varying> varyings;
	for (auto it = attributes.begin(); it != attributes.end(); it++) {
		varyings.push_back(vertexShader(uniform, *it));
	}
	return varyings;
}

vector<Triangle> assemblyStage(const vector<Varying> &varyings, const vector<unsigned int> &indices) {
	vector<Triangle> triangles;
	for (size_t i = 0; i < indices.size() - 2; i += 3) {
		Triangle triangle;
		triangle.a = varyings.at(indices[i]);
		triangle.b = varyings.at(indices[i + 1]);
		triangle.c = varyings.at(indices[i + 2]);
		// Viewport test, triangles outside [-1,-1],[1,1] are pruned
		if (
			all(greaterThanEqual(triangle.a.position.xy, vec2(-1.0f, -1.0f))) && all(lessThanEqual(triangle.a.position.xy, vec2(1.0f, 1.0f))) &&
			all(greaterThanEqual(triangle.b.position.xy, vec2(-1.0f, -1.0f))) && all(lessThanEqual(triangle.b.position.xy, vec2(1.0f, 1.0f))) &&
			all(greaterThanEqual(triangle.c.position.xy, vec2(-1.0f, -1.0f))) && all(lessThanEqual(triangle.c.position.xy, vec2(1.0f, 1.0f)))
		) {
			triangles.push_back(triangle);
		}
	}
	return triangles;
}

void fragmentStage(sampler2D &frameBuffer, const Uniform &uniform, const vector<Triangle> &triangles) {
	for (auto it = triangles.begin(); it != triangles.end(); it++) {
		// Prepare the input values
		const Triangle &triangle = *it;
		vec3 vertexA = vec3(triangle.a.position.xy, 1.0f);
		vec3 vertexB = vec3(triangle.b.position.xy, 1.0f);
		vec3 vertexC = vec3(triangle.c.position.xy, 1.0f);
		// Ascertain that no colinear vectors are encountered
		if (vertexA != vertexB && vertexA != vertexC && vertexB != vertexC) {
			mat3 pos2fb = scalingMatrix(vec3(vec2(frameBuffer.size), 1.0f)) * scalingMatrix(0.5f, 0.5f, 1.0f) * translationMatrix(1.0f, 1.0f) * scalingMatrix(1.0f, -1.0f, 1.0f);
			// Determine the fragment coordinates for all three vertices
			ivec2 fragCoordA = ivec2(vec2(pos2fb * vertexA));
			ivec2 fragCoordB = ivec2(vec2(pos2fb * vertexB));
			ivec2 fragCoordC = ivec2(vec2(pos2fb * vertexC));
			// Calculate the bounding box of the triangle
			ivec2 topLeft = min(min(fragCoordA, fragCoordB), fragCoordC);
			ivec2 bottomRight = max(max(fragCoordA, fragCoordB), fragCoordC);
			// Calculate the triangle vector space transform, relative to A
			mat2 invSpace = inverse(mat2(vec2(fragCoordB - fragCoordA), vec2(fragCoordC - fragCoordA)));
			//cout << "Scanning " << topLeft << "-" << bottomRight << endl;
			for (int y = topLeft.y; y < bottomRight.y; y++) {
				for (int x = topLeft.x; x < bottomRight.x; x++) {
					ivec2 fragPoint = ivec2(x, y);
					// Transform screen space coordinates into triangle space coordinates, relative to A
					vec2 trianglePoint = invSpace * vec2(fragPoint - fragCoordA);
					// Determine if the point is inside unitary triangle space (i.e. inside the triangle)
					if (trianglePoint.x >= 0.0f && trianglePoint.y >= 0.0f && trianglePoint.x + trianglePoint.y <= 1.0f) {
						// Calculate the squared relative distance to each vertex
						vec2 pointA = vec2(fragCoordA - fragPoint);
						vec2 pointB = vec2(fragCoordB - fragPoint);
						vec2 pointC = vec2(fragCoordC - fragPoint);
						//vec3 distance = vec3(dot(pointA, pointA), dot(pointB, pointB), dot(pointC, pointC));
						vec3 distance = vec3(length(pointA), length(pointB), length(pointC));
						// Calculate the varying blend from the distances
						float distanceSum = distance.x + distance.y + distance.z;
						vec3 bias = distance / distanceSum;
						//cout << " bias=" << bias << endl;
						Varying varying = blend(triangle.a, triangle.b, triangle.c, bias);
						// Execute the fragment shader for the varying blend
						Output fragment = fragmentShader(uniform, varying);
						// Write the result to the framebuffer
						texelSet(frameBuffer, fragPoint, 0, fragment.fragColor);
					}
				}
			}
		}
	}
}

void writeSamplerPnm(const sampler2D &sampler, const char *file) {
	ofstream out(file);
	if (!!out) {
		out << "P6" << std::endl;
		out << sampler.size.x << " " << sampler.size.y << std::endl;
		out << "255" << std::endl;
		for (unsigned int i = 0; i < sampler.size.x * sampler.size.y; i++) {
			vec4 rgba = sampler.samples[i];
			Vector<char, 4> bytes(ivec4(rgba * vec4(255.0f)));
			out.write(bytes.internal(), 3);
		}
		out.close();
	}
}

int main(int argc, char **argv) {
	Uniform uniform;
	uniform.Kd = vec4(1.0f);
	uniform.Ks = vec4(1.0f, 1.0f, 1.0f, 60.0f);
	uniform.Ka = vec3(0.2f);
	uniform.lightColor = vec3(1.0f, 1.0f, 0.8f);
	uniform.ambientColor = vec3(1.0f, 1.0f, 0.8f);;
	uniform.texture_d = header2sampler(sky_d_width, sky_d_height, sky_d_samples);
	uniform.texture_n = header2sampler(sky_n_width, sky_n_height, sky_n_samples);
	
	sampler2D frameBuffer;
	frameBuffer.size = ivec2(256, 256);
	frameBuffer.samples = new vec4[frameBuffer.size.x * frameBuffer.size.y];
	
	mat4 projectionMatrix = perspectiveMatrix(45.0f, float(frameBuffer.size.y) / float(frameBuffer.size.x), 0.1f, 10.0f);
	//mat4 objectMatrix = translationMatrix(0.0f, 0.0f, -2.0f);
	mat4 objectMatrix = rotationMatrix(radians(90.0f), 1.0f, 0.0f, 0.0f);
	vec3 lightPosition = vec3(10.0f, 10.0f, -20.0f);
	
	uniform.projectionModelViewMatrix = projectionMatrix * objectMatrix;
	uniform.normalMatrix = mat4(1.0f);
	uniform.objectLightPosition = vec3(inverse(objectMatrix) * vec4(lightPosition, 0.0f));
	
	vector<unsigned int> tpatches(teapot_patches, &teapot_patches[sizeof(teapot_patches) / sizeof(teapot_patches[0])]);
	vector<float> tvertices(teapot_vertices, &teapot_vertices[sizeof(teapot_vertices) / sizeof(teapot_vertices[0])]);
	vector<unsigned int> indices;
	vector<Attribute> attributes;
	bezier2mesh(tpatches, tvertices, 0, indices, attributes);
	cout << "indices=";
	bool first = true;
	for (auto it = indices.begin(); it != indices.end(); it++) {
		if (first) {
			first = false;
		} else {
			cout << ',';
		}
		cout << *it;
	}
	cout << endl;
	cout << "attributes=";
	first = true;
	for (auto it = attributes.begin(); it != attributes.end(); it++) {
		if (first) {
			first = false;
		} else {
			cout << ',';
		}
		cout << '(' << it->position << ',' << it->texc << ',' << it->normal << ',' << it->tangent << ')';
	}
	cout << endl;
	
	vector<Varying> varyings = vertexStage(uniform, attributes);
	vector<Triangle> triangles = assemblyStage(varyings, indices);
	fragmentStage(frameBuffer, uniform, triangles);
	
	if (argc >= 2) {
		writeSamplerPnm(frameBuffer, argv[1]);
	}
	
	return 0;
}
