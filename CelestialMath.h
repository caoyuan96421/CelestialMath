/*
 * CelestialMath.h
 *
 *  Created on: Feb 21, 2018
 *      Author: Yuan
 */

#ifndef CELESTIALMATH_H_
#define CELESTIALMATH_H_

#include <time.h>

struct EquatorialCoordinates {
	double dec;		// Declination
	double ra;		// Right ascension
	EquatorialCoordinates(double d = 0, double r = 0) :
			dec(d), ra(r) {
	}
};

struct LocalEquatorialCoordinates {
	double dec;		// Declination
	double ha;		// Hour angle
	LocalEquatorialCoordinates(double d = 0, double h = 0) :
			dec(d), ha(h) {
	}
};

struct AzimuthalCoordinates {
	double alt;		// Altitude
	double azi;		// Azimuth
	AzimuthalCoordinates(double a1 = 0, double a2 = 0) :
			alt(a1), azi(a2) {
	}
};

struct LocationCoordinates {
	double lat;		// Latitude
	double lon;		// Longtitude
	LocationCoordinates(double l1 = 0, double l2 = 0) :
			lat(l1), lon(l2) {
	}
};

/**
 * Utility functions for doing math on coordinates of the celestial sphere
 */
class CelestialMath {
public:
	CelestialMath() {
	}
	~CelestialMath() {
	}

	/*Basic conversion between reference frames*/
	static AzimuthalCoordinates localEquatorialToAzimuthal(LocalEquatorialCoordinates a, LocationCoordinates loc);
	static LocalEquatorialCoordinates azimuthalToLocalEquatorial(AzimuthalCoordinates b, LocationCoordinates loc);
	static double getGreenwichMeanSiderealTime(time_t timestamp);
	static double getLocalSiderealTime(time_t timestamp, LocationCoordinates loc);
	static LocalEquatorialCoordinates equatorialToLocalEquatorial(EquatorialCoordinates e, time_t timestamp, LocationCoordinates loc);
	static EquatorialCoordinates localEquatorialToEquatorial(LocalEquatorialCoordinates a, time_t timestamp, LocationCoordinates loc);

	/*Misalignment correction functions*/
	static LocalEquatorialCoordinates misalignedPolarAxis(LocalEquatorialCoordinates a, AzimuthalCoordinates mpa, LocationCoordinates loc);
};

#endif /* CELESTIALMATH_H_ */
