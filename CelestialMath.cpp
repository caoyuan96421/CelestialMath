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
#define M_PI 3.141592653589793
#endif

static const double RADIAN = 180.0 / M_PI;
static const double DEGREE = M_PI / 180.0;
static const double sidereal_day = 86164.09;

AzimuthalCoordinates CelestialMath::localEquatorialToAzimuthal(LocalEquatorialCoordinates a, LocationCoordinates loc) {
	//              cphi             lambda             lambda0
	double c1 = cos(a.ha * DEGREE), c2 = cos(a.dec * DEGREE), c3 = cos(loc.lat * DEGREE);
	double s1 = sin(a.ha * DEGREE), s2 = sin(a.dec * DEGREE), s3 = sin(loc.lat * DEGREE);
	double s4 = c1 * c2 * c3 + s2 * s3;
	s4 = clamp(s4);
	double y1 = s1 * c2, x1 = s2 * c3 - c1 * c2 * s3;

	AzimuthalCoordinates b;
	b.alt = asin(s4) * RADIAN;
	b.azi = atan2(y1, x1) * RADIAN;
	return b;
}

LocalEquatorialCoordinates CelestialMath::azimuthalToLocalEquatorial(AzimuthalCoordinates b, LocationCoordinates loc) {
	//              mu             eps             lambda0
	double c1 = cos(b.azi * DEGREE), c2 = cos(b.alt * DEGREE), c3 = cos(loc.lat * DEGREE);
	double s1 = sin(b.azi * DEGREE), s2 = sin(b.alt * DEGREE), s3 = sin(loc.lat * DEGREE);
	double s4 = c1 * c2 * c3 + s2 * s3;
	s4 = clamp(s4);
	double y1 = s1 * c2, x1 = s2 * c3 - c1 * c2 * s3;

	LocalEquatorialCoordinates a;
	a.dec = asin(s4) * RADIAN;
	a.ha = atan2(y1, x1) * RADIAN;
	return a;
}

double CelestialMath::getGreenwichMeanSiderealTime(time_t timestamp) {
	double jd = (double) timestamp * 1.1574074074074E-5 + 2440587.5; // Julian Date (J2000)
	double gmst = 280.46061837 + 360.985647366 * jd; // Greenwich mean sidereal time (angle)
	return remainder(gmst, 360.0);
}

double CelestialMath::getLocalSiderealTime(time_t timestamp, LocationCoordinates loc) {
	double gmst = getGreenwichMeanSiderealTime(timestamp);
	double lst = gmst + loc.lon * 1.00273790935; // Local sidereal time (angle)
	return remainder(lst, 360.0);
}

LocalEquatorialCoordinates CelestialMath::equatorialToLocalEquatorial(EquatorialCoordinates e, time_t timestamp, LocationCoordinates loc) {
	// From phi to cphi
	return EquatorialCoordinates(e.dec, remainder(getLocalSiderealTime(timestamp, loc) - e.ra, 360.0));
}

EquatorialCoordinates CelestialMath::localEquatorialToEquatorial(LocalEquatorialCoordinates a, time_t timestamp, LocationCoordinates loc) {
	// From cphi to phi
	return LocalEquatorialCoordinates(a.dec, remainder(getLocalSiderealTime(timestamp, loc) - a.ha, 360.0));
}

LocalEquatorialCoordinates CelestialMath::misalignedPolarAxis(LocalEquatorialCoordinates a, AzimuthalCoordinates mpa, LocationCoordinates loc) {
}
