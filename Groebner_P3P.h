// Copyright (C) 2018 Atsuhiko Banno (AIST)
//	atsuhiko.banno@aist.go.jp

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// If you use this program, please refer to the following paper:
// Atsuhiko Banno (AIST)
// "A P3P problem solver representing all parameters as a linear combination"
// Image Vision and Computing, 70 (2018) 55-62.


#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <Eigen/SVD>

class CGroebner
{
public:
	CGroebner(void);
	~CGroebner(void);
		
	int CalculatePosition(Eigen::Matrix3d globalPositions, Eigen::Matrix3d localViews, std::vector<Eigen::Matrix3d> &camR, std::vector<Eigen::Vector3d> &camT);
	// input 
	// globalPositions = | X1  X2  X3 |
	//                   | Y1  Y2  Y3 |
	//                   | Z1  Z2  Z3 |
	// localViews = | x1  x2  x3 |
	//              | y1  y2  y3 |  
	//              | z1  z2  z3 |

	// output
	// at most four sets of (R, T)

	// the notations are as in the above reference paper.


private:
	int solveQuartic(Eigen::VectorXd factors, Eigen::VectorXd & realRoots);
	int solveCubic(Eigen::VectorXd factors, Eigen::VectorXd & realRoots);

	inline double SQR(double x) {
		return x*x;
	};
};

