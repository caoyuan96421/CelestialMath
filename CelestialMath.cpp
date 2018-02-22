/*
 * CelestialMath.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: Yuan
 */

#include "CelestialMath.h"
#include <math.h>
#include <stdio.h>

static inline double clamp(double x) {
	return (x > 1) ? 1 : ((x < -1) ? -1 : x);
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double RADIAN = 180.0 / M_PI;
static const double DEGREE = M_PI / 180.0;
static const double sidereal_day = 86164.09;

CartesianVector CartesianVector::operator*(const Transformation &t) {
	return CartesianVector(t.a11 * x + t.a21 * y + t.a31 * z, t.a12 * x + t.a22 * y + t.a32 * z, t.a13 * x + t.a23 * y + t.a33 * z);
}

CartesianVector Transformation::operator*(const CartesianVector &vec) { // Left-product of matrix and vector
	return CartesianVector(a11 * vec.x + a12 * vec.y + a13 * vec.z, a21 * vec.x + a22 * vec.y + a23 * vec.z, a31 * vec.x + a32 * vec.y + a33 * vec.z);
}

AzimuthalCoordinates CelestialMath::localEquatorialToAzimuthal(const LocalEquatorialCoordinates &a, const LocationCoordinates &loc) {
	//              cphi             lambda             lambda0
	double c1 = cos(a.ha * DEGREE), c2 = cos(a.dec * DEGREE), c3 = cos(loc.lat * DEGREE);
	double s1 = sin(a.ha * DEGREE), s2 = sin(a.dec * DEGREE), s3 = sin(loc.lat * DEGREE);
	double s4 = c1 * c2 * c3 + s2 * s3;
	s4 = clamp(s4);
	double y1 = s1 * c2, x1 = s2 * c3 - c1 * c2 * s3;

	return AzimuthalCoordinates(asin(clamp(s4)) * RADIAN, atan2(y1, x1) * RADIAN);
}

LocalEquatorialCoordinates CelestialMath::azimuthalToLocalEquatorial(const AzimuthalCoordinates &b, const LocationCoordinates &loc) {
	//              mu             eps             lambda0
	double c1 = cos(b.azi * DEGREE), c2 = cos(b.alt * DEGREE), c3 = cos(loc.lat * DEGREE);
	double s1 = sin(b.azi * DEGREE), s2 = sin(b.alt * DEGREE), s3 = sin(loc.lat * DEGREE);
	double s4 = c1 * c2 * c3 + s2 * s3;
	double y1 = s1 * c2, x1 = s2 * c3 - c1 * c2 * s3;

	return LocalEquatorialCoordinates(asin(clamp(s4)) * RADIAN, atan2(y1, x1) * RADIAN);
}

double CelestialMath::getGreenwichMeanSiderealTime(time_t timestamp) {
	double jd = (double) timestamp * 1.1574074074074E-5 + 2440587.5; // Julian Date (J2000)
	double gmst = 280.46061837 + 360.985647366 * jd; // Greenwich mean sidereal time (angle)
	return remainder(gmst, 360.0);
}

double CelestialMath::getLocalSiderealTime(time_t timestamp, const LocationCoordinates &loc) {
	double gmst = getGreenwichMeanSiderealTime(timestamp);
	double lst = gmst + loc.lon * 1.00273790935; // Local sidereal time (angle)
	return remainder(lst, 360.0);
}

LocalEquatorialCoordinates CelestialMath::equatorialToLocalEquatorial(const EquatorialCoordinates &e, time_t timestamp, const LocationCoordinates &loc) {
	// From phi to cphi
	return LocalEquatorialCoordinates(e.dec, remainder(getLocalSiderealTime(timestamp, loc) - e.ra, 360.0));
}

EquatorialCoordinates CelestialMath::localEquatorialToEquatorial(const LocalEquatorialCoordinates &a, time_t timestamp, const LocationCoordinates &loc) {
	// From cphi to phi
	return EquatorialCoordinates(a.dec, remainder(getLocalSiderealTime(timestamp, loc) - a.ha, 360.0));
}

void CelestialMath::getMisalignedPolarAxisTransformation(Transformation *t, const AzimuthalCoordinates &p, const LocationCoordinates &loc) {
	double c1 = cos(p.azi * DEGREE), c2 = cos(p.alt * DEGREE), c3 = cos(loc.lat * DEGREE);
	double s1 = sin(p.azi * DEGREE), s2 = sin(p.alt * DEGREE), s3 = sin(loc.lat * DEGREE);
	// Matrix to convert from basis vectors in misaligned PA to correct PA
	t->a11 = c1 * s2 * s3 + c2 * c3;
	t->a12 = -s1 * s3;
	t->a13 = -c1 * c2 * s3 + s2 * c3;
	t->a21 = s1 * s2;
	t->a22 = c1;
	t->a23 = -s1 * c2;
	t->a31 = -c1 * s2 * c3 + c2 * s3;
	t->a32 = s1 * c3;
	t->a33 = c1 * c2 * c3 + s2 * s3;
}

LocalEquatorialCoordinates CelestialMath::transform(Transformation* t, const LocalEquatorialCoordinates& a) {
	double c1 = cos(a.dec * DEGREE), c2 = cos(a.ha * DEGREE);
	double s1 = sin(a.dec * DEGREE), s2 = cos(a.ha * DEGREE);
	CartesianVector X = CartesianVector(c1 * c2, -c1 * s2, s1) * (*t);

	return LocalEquatorialCoordinates(asin(clamp(X.z)) * RADIAN, atan2(-X.y, X.x) * RADIAN);
}
