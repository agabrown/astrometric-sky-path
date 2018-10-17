/*
 * astrometry.h
 *
 * Include file that defines convenient data structures.
 *
 * Anthony Brown 23.08.2011 - Oct 2018
 */
typedef struct {

/*
 * Define a data structure describing a single star in terms of its astrometry and corresponding phase
 * space coordinates.
 */
  double rightAscension; /* Right ascension in radians */
  double declination;    /* Declination in radians */
  double parallax;        /* Parallax (mas) */
  double properMotionRa; /* Proper motion in right ascension (mu_alpha*cos(delta), mas/yr) */
  double properMotionDec; /* Proper motion in declination (mas/yr) */
  double radialVelocity;  /* Radial Velocity (km/s) */
  double phaseSpaceCoordinates[2][3]; /* Phase space coordinates [0][0-2] = (x,y,z) (pc) */
                                      /*                         [1][0-2] = (vx,vy,vz) (km/s) */
}
star;

/*
 * Data structure that holds the observation time interval and the reference epoch.
 */
typedef struct{
    double refEpJ; // Reference epoch for the astrometric parameters (Julian year).
    double beginEpJ; // Start of observation itme interval (Julian year).
    double deltaEpJ; // Length of observation time interval in Julian years.
}
obsEpochs;
