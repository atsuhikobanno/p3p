// Copyright (C) 2018 Atsuhiko Banno (AIST)

// Groebner_P3P is free software : you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Groebner_P3P is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Groebner_P3P. If not, see <http://www.gnu.org/licenses/>.

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

