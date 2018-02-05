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


#include "Groebner_P3P.h"

using namespace Eigen;

CGroebner::CGroebner(void)
{
}

CGroebner::~CGroebner(void)
{
}

int CGroebner::CalculatePosition(Matrix3d globalPositions, Matrix3d localViews, std::vector<Matrix3d> &camR, std::vector<Vector3d> &camT)
{
//	Coordinate changes
	double BaseScale = 10.0;
	Matrix3d R, Rg, Rc, tM;
	Vector3d tmpV, nx, ny, nz;

	globalPositions = (1 / BaseScale) * globalPositions;

	nx = globalPositions.col(1) - globalPositions.col(0);
	nx.normalize();

	tmpV = globalPositions.col(2) - globalPositions.col(0);
	nz = nx.cross(tmpV);
	nz.normalize();

	ny = nz.cross(nx);

	Rg.col(0) = nx;
	Rg.col(1) = ny;
	Rg.col(2) = nz;

	tM.col(0) = Vector3d::Zero();
	tM.col(1) = globalPositions.col(1) - globalPositions.col(0);
	tM.col(2) = globalPositions.col(2) - globalPositions.col(0);

	tM = Rg.transpose() * tM;

	std::vector<Vector3d> neo_globalPositions(3);
	neo_globalPositions[0] << tM(0, 0), tM(1, 0), tM(2, 0);
	neo_globalPositions[1] << tM(0, 1), tM(1, 1), tM(2, 1);
	neo_globalPositions[2] << tM(0, 2), tM(1, 2), tM(2, 2);

	VectorXd r0 = VectorXd::Zero(12);
	VectorXd r1 = VectorXd::Zero(12);
	VectorXd r2 = VectorXd::Zero(12);

	r0[ 9] = 1;
	r1[10] = 1;
	r2[11] = 1;

	double XY = 1 / neo_globalPositions[1].x() / neo_globalPositions[2].y();
	double XmY = (neo_globalPositions[1].x() - neo_globalPositions[2].x()) * XY;

	localViews.col(0).normalize();
	localViews.col(1).normalize();
	localViews.col(2).normalize();

	r0[0] = - localViews(0, 0) / neo_globalPositions[1].x();
	r0[1] = - localViews(0, 0) * XmY;
	r0[2] = - localViews(1, 0) / neo_globalPositions[1].x();
	r0[3] = - localViews(1, 0) * XmY;
	r0[4] = - localViews(2, 0) / neo_globalPositions[1].x();
	r0[5] = - localViews(2, 0) * XmY;

	r1[0] =   localViews(0, 1) / neo_globalPositions[1].x();
	r1[1] = - localViews(0, 1) * neo_globalPositions[2].x() * XY;
	r1[2] =   localViews(1, 1) / neo_globalPositions[1].x();
	r1[3] = - localViews(1, 1) * neo_globalPositions[2].x() * XY;
	r1[4] =   localViews(2, 1) / neo_globalPositions[1].x();
	r1[5] = - localViews(2, 1) * neo_globalPositions[2].x() * XY;

	r2[1] = localViews(0, 2) / neo_globalPositions[2].y();
	r2[3] = localViews(1, 2) / neo_globalPositions[2].y();
	r2[5] = localViews(2, 2) / neo_globalPositions[2].y();

	r0[6] = localViews(0, 0);
	r0[7] = localViews(1, 0);
	r0[8] = localViews(2, 0);

	double c12 = localViews.col(0).dot(localViews.col(1));
	double c23 = localViews.col(1).dot(localViews.col(2));
	double c31 = localViews.col(2).dot(localViews.col(0));

	double y_32 = SQR(neo_globalPositions[2].y());
	double r_12 = neo_globalPositions[1].x() - neo_globalPositions[2].x();
	double f1 = y_32 - r_12 * r_12;
	double f2 = -2.0*(y_32 + r_12*neo_globalPositions[2].x()) * c12;
	double f3 = 2.0*neo_globalPositions[1].x() * r_12*c31;
	double f4 = y_32 - neo_globalPositions[2].x() * neo_globalPositions[2].x();
	double f5 = 2.0*neo_globalPositions[1].x() * neo_globalPositions[2].x() * c23;
	double f6 = - neo_globalPositions[1].x() * neo_globalPositions[1].x();

	double g1 = r_12;
	double g2 = (neo_globalPositions[2].x() - r_12)*c12;
	double g3 = -neo_globalPositions[1].x() * c31;
	double g4 = -neo_globalPositions[2].x();
	double g5 = neo_globalPositions[1].x() * c23;

	double f11 = f1*f1;
	double f12 = f1*f2;
	double f13 = f1*f3;
	double f14 = f1*f4;
	double f15 = f1*f5;
	double f16 = f1*f6;

	double f22 = f2*f2;
	double f23 = f2*f3;
	double f24 = f2*f4;
	double f25 = f2*f5;
	double f26 = f2*f6;

	double f33 = f3*f3;
	double f34 = f3*f4;
	double f35 = f3*f5;
	double f36 = f3*f6;

	double f44 = f4*f4;
	double f45 = f4*f5;
	double f46 = f4*f6;

	double f55 = f5*f5;
	double f56 = f5*f6;

	double f66 = f6*f6;

	double g11 = g1*g1;
	double g12 = g1*g2;
	double g13 = g1*g3;
	double g14 = g1*g4;
	double g15 = g1*g5;
		;
	double g22 = g2*g2;
	double g23 = g2*g3;
	double g24 = g2*g4;
	double g25 = g2*g5;

	double g33 = g3*g3;
	double g34 = g3*g4;
	double g35 = g3*g5;

	double g44 = g4*g4;
	double g45 = g4*g5;

	double g55 = g5*g5;

// 4th-order equation 
	int sol=0;
	VectorXd realRoots(4);
	VectorXd factors(5);

	factors[0] = f44*g11 + f11*g44 - f24*g12 + f22*g14 + f14*g22 - f12*g24 - 2*f14*g14;

	factors[1] = 2*(f45*g11 + f11*g45 + f23*g14 + f14*g23 - f15*g14 - f14*g15) 
				- f34*g12 - f25*g12 - f24*g13 + f22*g15 + f15*g22 - f13*g24 - f12*g34 - f12*g25;

	factors[2] = 2*(f46*g11 + f23*g15 + f15*g23 - f16*g14 - f15*g15) + f55*g11 + f11*g55 - f35*g12
				- f26*g12 - f34*g13 - f25*g13 + f33*g14 + f16*g22 + f14*g33 - f13*g34 - f13*g25 - f12*g35;
	factors[3] = 2*(f56*g11 + f16*g23 - f16*g15) - f36*g12 - f35*g13 - f26*g13 + f33*g15 + f15*g33 - f13*g35;
	factors[4] = f66*g11 - f36*g13 + f16*g33;

	if(fabs(factors[0])<1.0e-12)
		sol = this->solveCubic(factors, realRoots);
	else
		sol = this->solveQuartic(factors, realRoots);

	if(sol==0)
		return false;

	int available=0;

	double x, y, z, nrm;
	Vector3d tc, T;

	camR.clear();
	camT.clear();
	camR.reserve(sol);
	camT.reserve(sol);

	for(int i=0; i<sol; i++){
		y = realRoots[i];

		x = f6*g1 + y*( f5*g1 - f1*g5 + y*(f4*g1 - f1*g4));
		x /= (f1*g2 - f2*g1)*y + f1*g3 - f3*g1;

		if(x<0 || y<0)
			continue;

		nrm = SQR(x*r0[0] + y*r1[0]) + SQR(x*r0[2] + y*r1[2]) + SQR(x*r0[4] + y*r1[4])
			+ SQR(x*r0[1] + y*r1[1] + r2[1]) + SQR(x*r0[3] + y*r1[3] + r2[3]) + SQR(x*r0[5] + y*r1[5] + r2[5]);
		z = 1/sqrt(nrm/2);
		x *= z;
		y *= z;

		Vector3d R1, R2, R3;
		R1 << x*r0[0] + y*r1[0], x*r0[2] + y*r1[2], x*r0[4] + y*r1[4];
		R2 << x*r0[1] + y*r1[1] + z*r2[1], x*r0[3] + y*r1[3] + z*r2[3], x*r0[5] + y*r1[5] + z*r2[5];
		R1.normalize();
		R2.normalize();
		R3 << R1.cross(R2);
		R3.normalize();

		Rc.row(0) = R1.transpose();
		Rc.row(1) = R2.transpose();
		Rc.row(2) = R3.transpose();

		//	Modify as a rotation matrix
		JacobiSVD< Matrix3d > svd(Rc, ComputeFullU | ComputeFullV);
		double det = (svd.matrixU() * svd.matrixV().transpose()).determinant();
		Matrix3d Id = Matrix3d::Identity();
		Id(2, 2) = det;
		Rc = svd.matrixU() * Id * svd.matrixV().transpose();

		tc << x*r0[6], x*r0[7], x*r0[8];

		R = Rg * Rc;
		T = globalPositions.col(0) - R * tc;
		T = BaseScale * T;

		if(R.trace()+1<=0)
			continue;

		camR.push_back(R);
		camT.push_back(T);

		available++;
	}
	
	return available;
}

int CGroebner::solveQuartic(VectorXd factors, VectorXd & realRoots)
{
	double A = factors[0];
	double B = factors[1];
	double C = factors[2];
	double D = factors[3];
	double E = factors[4];

	double A_pw2 = A*A;
	double B_pw2 = B*B;
	double A_pw3 = A_pw2*A;
	double B_pw3 = B_pw2*B;
	double A_pw4 = A_pw3*A;
	double B_pw4 = B_pw3*B;

	double alpha = -3*B_pw2/(8*A_pw2)+C/A;
	double beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
	double gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

	double alpha_pw2 = alpha*alpha;
	double alpha_pw3 = alpha_pw2*alpha;

	std::complex<double> P (-alpha_pw2/12 - gamma,0);
	std::complex<double> Q (-alpha_pw3/108 + alpha*gamma/3 - beta*beta/8,0);
	std::complex<double> R = -Q/2.0+sqrt(Q*Q/4.0 + P*P*P/27.0);

	std::complex<double> U = pow(R,(1.0/3.0));
	std::complex<double> y;

	if (U.real() == 0)
		y = -5.0*alpha/6.0-pow(Q,(1.0/3.0));
	else
		y = -5.0*alpha/6.0-P/(3.0*U)+U;

	std::complex<double> w = sqrt(alpha+2.0*y);

	std::complex<double> temp;

	int sol = 0;
	temp = -B/(4.0*A) + 0.5*(w+sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
	realRoots[sol++] = temp.real();
	temp = -B/(4.0*A) + 0.5*(w-sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
	realRoots[sol++] = temp.real();
	temp = -B/(4.0*A) + 0.5*(-w+sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
	realRoots[sol++] = temp.real();
	temp = -B/(4.0*A) + 0.5*(-w-sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
	realRoots[sol++] = temp.real();

	return sol;
}

int CGroebner::solveCubic(VectorXd factors, VectorXd &realRoots)
{
	double A = factors[1];
	double B = factors[2];
	double C = factors[3];
	double D = factors[4];

	std::complex<double> w(-0.5, sqrt(3)/2);
	std::complex<double> w2 = w*w;

	double alpha, beta;

	alpha = 3*(27*A*A*D*D - 18*A*B*C*D + 4*A*C*C*C + 4*B*B*B*D - B*B*C*C);
	alpha = 3*A*sqrt(alpha);

	beta = -27*A*A*D + 9*A*B*C - 2*B*B*B; 

	std::complex<double> P, Q, temp;
	
	int sol = 0;

	P = std::complex<double>(1, 0);
	Q = std::complex<double>(1, 0);
	temp = (-2*B + P*pow(4*(beta + alpha), 1/3) + Q*pow(4*(beta - alpha), 1/3))/(6*A);
	realRoots[sol++] = temp.real();

	P = w;
	Q = w2;
	temp = (-2*B + P*pow(4*(beta + alpha), 1/3) + Q*pow(4*(beta - alpha), 1/3))/(6*A);
	realRoots[sol++] = temp.real();

	P = w2;
	Q = w;
	temp = (-2*B + P*pow(4*(beta + alpha), 1/3) + Q*pow(4*(beta - alpha), 1/3))/(6*A);
	realRoots[sol++] = temp.real();

	return sol;
}
