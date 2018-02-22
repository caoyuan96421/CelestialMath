/*
 * main.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: Yuan
 */

#include <stdio.h>
#include "CelestialMath.h"
#include <math.h>

int test1(double dec, double ha) {
	LocalEquatorialCoordinates coord(dec, ha);
	LocationCoordinates loc(42.0, 0);
	AzimuthalCoordinates ac = CelestialMath::localEquatorialToAzimuthal(coord,
			loc);
	printf("Alt = %f,\tAzi = %f\n", ac.alt, ac.azi);
	coord = CelestialMath::azimuthalToLocalEquatorial(ac, loc);
	printf("Dec = %.10f,\tHA = %.10f\n", coord.dec, coord.ha);

	return (fabs(coord.dec - dec) < 1e-7
			&& (fabs(dec) > 90.0 - 1e-7
					|| fabs(remainder(coord.ha - ha, 360.0)) < 1e-7));
}

int main() {
	printf("%d\n", test1(0, 0));
	printf("%d\n", test1(90, 0));
	printf("%d\n", test1(-90, 0));
	printf("%d\n", test1(35, 360));
	printf("%d\n", test1(46.12345, 720));
	printf("%d\n", test1(89.99998, 187));
	return 0;
}

