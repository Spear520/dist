#pragma once
#include <cmath>
#include <algorithm>
#include <iostream>
#include"point_line.h"


// point-triangle derivatives; returns squared distance
// input: x[12] is the array of coordinates of 4 points p0, p1, p2, p3
// p0 is the point from which the distance is computed,
// p1,p2,p3 form the triangle
// the function returns the squared distance
// fd is the output array of first derivatives of the squared distance
// sd is the output array of the second derivatives of the squared distance
// zeta2 and zeta3 are barycentric coordinates of the closest point on the triangle
// zeta1 can be obtained as 1-(zeta2+zeta3)
double pt(double(&x)[12], double(&fd)[12], double(&sd)[12][12], double &zeta2, double &zeta3);

// POINT-PLANE (intended to use in the interior of the triangle)
// second derivatives of a
double a2[12][12] = {
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0 },
	{ 0, 0, 0, -2, 0, 0, 2, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };

// second derivatives of b
double b2[12][12] = {
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 2, 0, 0, -1, 0, 0, -1, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, -1, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 0, -1, 0, 0, -1 },
	{ 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0 },
	{ 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0 },
	{ 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1 },
	{ 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0 } };

// second derivatives of c
double c2[12][12] = {
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2, 0 },
	{ 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, -2 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0 },
	{ 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0 },
	{ 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 2 } };

// second derivatives of d
double d2[12][12] = {
	{ 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0 },
	{ 1, 0, 0, -2, 0, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 1, 0, 0, -2, 0, 0, 1, 0, 0, 0, 0 },
	{ 0, 0, 1, 0, 0, -2, 0, 0, 1, 0, 0, 0 },
	{ -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };

// second derivatives of e
double e2[12][12] = {
	{ 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0 },
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0 },
	{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1 },
	{ 1, 0, 0, -2, 0, 0, 0, 0, 0, 1, 0, 0 },
	{ 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, 1, 0 },
	{ 0, 0, 1, 0, 0, -2, 0, 0, 0, 0, 0, 1 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0 } };

// second derivatives of f
double f2[12][12] = {
	{ 2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 0, 0 },
	{ -2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, -2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, -2, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };


// Kronecker delta
double xd(int idx1, int idx2) { return idx1 == idx2 ? 1. : 0; }

// input: x[12] array of p0,p1,p2,p3
// output: first derivatives of squared distance fd[12]
// second derivatives of squared distance sd[12][12]
// return value is the squared distance
// note: don't call this function directly
double point_plane_distance(double(&x)[12], double(&fd)[12], double(&sd)[12][12])
{
	double output_s, output_t; // for testing
	double x0 = x[0];
	double x1 = x[1];
	double x2 = x[2];
	double x3 = x[3];
	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];
	double x8 = x[8];
	double x9 = x[9];
	double x10 = x[10];
	double x11 = x[11];

	double abcdef[6] = { (-x3 + x6)*(-x3 + x6) + (-x4 + x7)*(-x4 + x7) + (-x5 + x8)*(-x5 + x8),
		(x10 - x4)*(-x4 + x7) + (x11 - x5)*(-x5 + x8) + (-x3 + x6)*(-x3 + x9),
		(x10 - x4)*(x10 - x4) + (x11 - x5)*(x11 - x5) + (-x3 + x9)*(-x3 + x9),
		(-x0 + x3)*(-x3 + x6) + (-x1 + x4)*(-x4 + x7) + (-x2 + x5)*(-x5 + x8),
		(x10 - x4)*(-x1 + x4) + (x11 - x5)*(-x2 + x5) + (-x0 + x3)*(-x3 + x9),
		(-x0 + x3)*(-x0 + x3) + (-x1 + x4)*(-x1 + x4) + (-x2 + x5)*(-x2 + x5) };

	double a = abcdef[0];
	double b = abcdef[1];
	double c = abcdef[2];
	double d = abcdef[3];
	double e = abcdef[4];
	double f = abcdef[5];

	double det = a*c - b*b;
	double detsq = det * det;
	double detcube = detsq * det;

	double s = b*e - c*d;
	double t = b*d - a*e;

	double invDet = 1. / det;
	s *= invDet;
	t *= invDet;
	output_s = s;
	output_t = t;
	double u = 1 - (s + t);

	double sqrDistance = (-x0 + x6*s + x9*t + x3*u)*(-x0 + x6*s + x9*t + x3*u) +
		(-x1 + x7*s + x10*t + x4*u)*(-x1 + x7*s + x10*t + x4*u) +
		(-x2 + x8*s + x11*t + x5*u)*(-x2 + x8*s + x11*t + x5*u);
	double dist = sqrt(sqrDistance);

	// select either normal case or degenerate case
	// u = zeta1; s = zeta2; t = zeta3;
	double s2[12][12], t2[12][12], det2[12][12];
	double u1[12], u2[12][12];

	// derivatives of s and t
	// first derivatives of the above quantities
	double a1[12] = { 0,0,0,-2 * (-x3 + x6),-2 * (-x4 + x7),-2 * (-x5 + x8),2 * (-x3 + x6),2 * (-x4 + x7),2 * (-x5 + x8),0,0,0 };
	double b1[12] = { 0,0,0,2 * x3 - x6 - x9,-x10 + 2 * x4 - x7,-x11 + 2 * x5 - x8,-x3 + x9,x10 - x4,x11 - x5,-x3 + x6,-x4 + x7,-x5 + x8 };
	double c1[12] = { 0,0,0,-2 * (-x3 + x9),-2 * (x10 - x4),-2 * (x11 - x5),0,0,0,2 * (-x3 + x9),2 * (x10 - x4),2 * (x11 - x5) };
	double d1[12] = { x3 - x6, x4 - x7, x5 - x8, x0 - 2 * x3 + x6, x1 - 2 * x4 + x7, x2 - 2 * x5 + x8, -x0 + x3, -x1 + x4, -x2 + x5, 0, 0, 0 };
	double e1[12] = { x3 - x9, -x10 + x4, -x11 + x5, x0 - 2 * x3 + x9, x1 + x10 - 2 * x4, x11 + x2 - 2 * x5, 0, 0, 0, -x0 + x3, -x1 + x4, -x2 + x5 };
	double f1[12] = { -2 * (-x0 + x3),-2 * (-x1 + x4),-2 * (-x2 + x5),2 * (-x0 + x3),2 * (-x1 + x4),2 * (-x2 + x5),0,0,0,0,0,0 };
	double s1[12], t1[12], det1[12];

	// first derivatives 
	for (int i = 0; i < 12; i++)
	{
		det1[i] = c*a1[i] + a*c1[i] - 2 * b*b1[i];
		s1[i] = ((c*d - b*e)*det1[i]) / detsq + ((e*b1[i] + b*e1[i]) - (d*c1[i] + c*d1[i])) / det;
		t1[i] = ((a*e - b*d)*det1[i]) / detsq + ((d*b1[i] + b*d1[i]) - (a*e1[i] + e*a1[i])) / det;
		u1[i] = -(s1[i] + t1[i]);

		fd[i] = -2 * (x0 - x6*s - x9*t - x3*u)*
			(x6*s1[i] + x9*t1[i] + x3*u1[i] - xd(0, i) + u*xd(3, i) +
				s*xd(6, i) + t*xd(9, i)) -
			2 * (x1 - x7*s - x10*t - x4*u)*
			(x7*s1[i] + x10*t1[i] + x4*u1[i] - xd(1, i) + u*xd(4, i) +
				s*xd(7, i) + t*xd(10, i)) -
			2 * (x2 - x8*s - x11*t - x5*u)*
			(x8*s1[i] + x11*t1[i] + x5*u1[i] - xd(2, i) + u*xd(5, i) +
				s*xd(8, i) + t*xd(11, i));
	}

	// loop may be simplified, because all matrices are symmetric
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
		{
			det2[i][j] = -2 * b1[i] * b1[j] + a1[j] * c1[i] + a1[i] * c1[j] + c*a2[i][j] - 2 * b*b2[i][j] + a*c2[i][j];

			s2[i][j] =
				+(-(c1[j] * d1[i]) - c1[i] * d1[j] + b1[j] * e1[i] + b1[i] * e1[j] + e*b2[i][j] - d*c2[i][j] - c*d2[i][j] + b*e2[i][j]) / det
				- ((det1[j] * (e*b1[i] - d*c1[i] - c*d1[i] + b*e1[i])) + (det1[i] * (e*b1[j] - d*c1[j] - c*d1[j] + b*e1[j])) + ((-(c*d) + b*e)*det2[i][j])) / detsq
				+ (2 * (-(c*d) + b*e)*det1[i] * det1[j]) / detcube;

			t2[i][j] =
				+(b1[j] * d1[i] + b1[i] * d1[j] - a1[j] * e1[i] - a1[i] * e1[j] - e*a2[i][j] + d*b2[i][j] + b*d2[i][j] - a*e2[i][j]) / det
				- ((det1[j] * (-(e*a1[i]) + d*b1[i] + b*d1[i] - a*e1[i])) + (det1[i] * (-(e*a1[j]) + d*b1[j] + b*d1[j] - a*e1[j])) + ((b*d - a*e)*det2[i][j])) / detsq
				+ (2 * (b*d - a*e)*det1[i] * det1[j]) / detcube;

			u2[i][j] = -(s2[i][j] + t2[i][j]);

			sd[i][j] = 2 * ((x6*s1[i] + x9*t1[i] + x3*u1[i] - xd(0, i) + u*xd(3, i) +
				s*xd(6, i) + t*xd(9, i))*
				(x6*s1[j] + x9*t1[j] + x3*u1[j] - xd(0, j) + u*xd(3, j) +
					s*xd(6, j) + t*xd(9, j)) -
					(x0 - x6*s - x9*t - x3*u)*
				(-0 + 0 * s + 0 * t + 0 * u + x6*s2[i][j] + x9*t2[i][j] +
					x3*u2[i][j] + u1[j] * xd(3, i) + u1[i] * xd(3, j) + s1[j] * xd(6, i) +
					s1[i] * xd(6, j) + t1[j] * xd(9, i) + t1[i] * xd(9, j)) +
					(x7*s1[i] + x10*t1[i] + x4*u1[i] - xd(1, i) + u*xd(4, i) +
						s*xd(7, i) + t*xd(10, i))*
						(x7*s1[j] + x10*t1[j] + x4*u1[j] - xd(1, j) + u*xd(4, j) +
							s*xd(7, j) + t*xd(10, j)) -
							(x1 - x7*s - x10*t - x4*u)*
				(-0 + 0 * s + 0 * t + 0 * u + x7*s2[i][j] + x10*t2[i][j] +
					x4*u2[i][j] + u1[j] * xd(4, i) + u1[i] * xd(4, j) + s1[j] * xd(7, i) +
					s1[i] * xd(7, j) + t1[j] * xd(10, i) + t1[i] * xd(10, j)) +
					(x8*s1[i] + x11*t1[i] + x5*u1[i] - xd(2, i) + u*xd(5, i) +
						s*xd(8, i) + t*xd(11, i))*
						(x8*s1[j] + x11*t1[j] + x5*u1[j] - xd(2, j) + u*xd(5, j) +
							s*xd(8, j) + t*xd(11, j)) -
							(x2 - x8*s - x11*t - x5*u)*
				(-0 + 0 * s + 0 * t + 0 * u + x8*s2[i][j] + x11*t2[i][j] +
					x5*u2[i][j] + u1[j] * xd(5, i) + u1[i] * xd(5, j) + s1[j] * xd(8, i) +
					s1[i] * xd(8, j) + t1[j] * xd(11, i) + t1[i] * xd(11, j)));
		}

	return sqrDistance;
}


// POINT-TRIANGE (ALL CASES COMBINED)
double pt(double(&x)[12], double(&fd)[12], double(&sd)[12][12], double &zeta2, double &zeta3)
{
	double x0 = x[0];
	double x1 = x[1];
	double x2 = x[2];
	double x3 = x[3];
	double x4 = x[4];
	double x5 = x[5];
	double x6 = x[6];
	double x7 = x[7];
	double x8 = x[8];
	double x9 = x[9];
	double x10 = x[10];
	double x11 = x[11];

	double a = (-x3 + x6)*(-x3 + x6) + (-x4 + x7)*(-x4 + x7) + (-x5 + x8)*(-x5 + x8);
	double b = (x10 - x4)*(-x4 + x7) + (x11 - x5)*(-x5 + x8) + (-x3 + x6)*(-x3 + x9);
	double c = (x10 - x4)*(x10 - x4) + (x11 - x5)*(x11 - x5) + (-x3 + x9)*(-x3 + x9);
	double d = (-x0 + x3)*(-x3 + x6) + (-x1 + x4)*(-x4 + x7) + (-x2 + x5)*(-x5 + x8);
	double e = (x10 - x4)*(-x1 + x4) + (x11 - x5)*(-x2 + x5) + (-x0 + x3)*(-x3 + x9);

	double det = a*c - b*b;
	double s = b*e - c*d;
	double t = b*d - a*e;
	double sqrDistance;

	if (s + t <= det)
	{
		if (s < 0)
		{
			if (t < 0)  // region 4
			{
				if (d < 0)
				{
					t = 0;
					if (-d >= a)
					{
						s = 1;
						sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
					}
					else {
						s = -d / a;
						sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 1, fd, sd);
					}
				}
				else {
					s = 0;
					if (e >= 0)
					{
						t = 0;
						sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
					}
					else if (-e >= c)
					{
						t = 1;
						sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
					}
					else {
						t = -e / c;
						sqrDistance = vertex_edge_distance_and_derivs_12(x, 3, 1, fd, sd);
					}
				}
			}
			else  // region 3
			{
				s = 0;
				if (e >= 0)
				{
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
				}
				else if (-e >= c)
				{
					t = 1;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
				}
				else {
					t = -e / c;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 3, 1, fd, sd);
				}
			}
		}
		else if (t < 0)  // region 5
		{
			t = 0;
			if (d >= 0)
			{
				s = 0;
				sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
			}
			else if (-d >= a)
			{
				s = 1;
				sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
			}
			else {
				s = -d / a;
				sqrDistance = vertex_edge_distance_and_derivs_12(x, 1, 2, fd, sd);
			}
		}
		else  // region 0
		{
			double invDet = (1) / det;
			s *= invDet;
			t *= invDet;
			sqrDistance = point_plane_distance(x, fd, sd);
		}
	}
	else {
		double tmp0, tmp1, numer, denom;

		if (s < 0)  // region 2
		{
			tmp0 = b + d;
			tmp1 = c + e;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2 * b + c;
				if (numer >= denom)
				{
					s = 1;
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
				}
				else {
					s = numer / denom;
					t = 1 - s;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 3, fd, sd);
				}
			}
			else {
				s = 0;
				if (tmp1 <= 0)
				{
					t = 1;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
				}
				else if (e >= 0)
				{
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
				}
				else {
					t = -e / c;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 1, 3, fd, sd);
				}
			}
		}
		else if (t < 0)  // region 6
		{
			tmp0 = b + e;
			tmp1 = a + d;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2 * b + c;
				if (numer >= denom)
				{
					t = 1;
					s = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
				}
				else {
					t = numer / denom;
					s = 1 - t;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 3, fd, sd);
				}
			}
			else {
				t = 0;
				if (tmp1 <= 0)
				{
					s = 1;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
				}
				else if (d >= 0)
				{
					s = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 1, x, fd, sd);
				}
				else {
					s = -d / a;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 1, 2, fd, sd);
				}
			}
		}
		else  // region 1
		{
			numer = c + e - b - d;
			if (numer <= 0)
			{
				s = 0;
				t = 1;
				sqrDistance = vertex_vertex_distance_and_derivs_12(0, 3, x, fd, sd);
			}
			else {
				denom = a - 2 * b + c;
				if (numer >= denom)
				{
					s = 1;
					t = 0;
					sqrDistance = vertex_vertex_distance_and_derivs_12(0, 2, x, fd, sd);
				}
				else {
					s = numer / denom;
					t = 1 - s;
					sqrDistance = vertex_edge_distance_and_derivs_12(x, 2, 3, fd, sd);
				}
			}
		}
	}

	zeta2 = s;
	zeta3 = t;
	return sqrDistance;
}