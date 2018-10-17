/*
 * Simplified calculation of the path of a star on the sky over time.
 *
 * Given the astrometric parameters (positions, parallax, proper motions, radial velocity) of a star
 * at a certain epoch, trace the path of the star on the sky as a function of time.
 *
 * Make use of the IAU SOFA library: http://www.iausofa.org/
 *
 * The sky path is obtained by transforming for each epoch the barycentric coordinate direction of
 * the star (centred on the Solar system barycentre) to topocentric coordinate direction (centred on
 * the Earth). This is done through the combination of the iauEpv00, iauPmpx, and iauC2s function
 * calls.
 *
 * The resulting coordinate directions do not contain the relativistic light deflection, abberation,
 * and other effects, and thus do not represent the actually observed proper directions for a
 * certain datei and place.
 *
 * The input parameters can be given either as phase space coordinates (x,y,z) and (vx,vy,vz) or as
 * the astrometric parameters plus radial velocity.
 *
 * Anthony Brown: 23 Aug 2011 (translated from the older Fortran version)
 *
 * Last updated: Oct 2018
 */

#include "skypath.h"

int main (int argc, char **argv) {

  double *times, *deltaRaCosDec, *deltaDec;
  star myStar, sofaStar;
  obsEpochs theObsEpochs;

  const int NPOINTS=1000;

  parseArgs(argc, argv, &myStar, &theObsEpochs);
  myUnitsToSofa(&myStar, &sofaStar);

  deltaRaCosDec = dvector(NPOINTS);
  deltaDec = dvector(NPOINTS);
  times = dvector(NPOINTS);
  calcSkyPath(&theObsEpochs, NPOINTS, &sofaStar, times, deltaRaCosDec, deltaDec);

  /*
   * Write out results.
   */
  writeToStdout(NPOINTS, times, deltaRaCosDec, deltaDec);

  free_dvector(deltaRaCosDec);
  free_dvector(deltaDec);
  free_dvector(times);

  return 0;
}

void writeToStdout(int npoints, double *times, double* deltaRaCosDec, double *deltaDec) {
  int k;
  for (k=0; k<npoints; k++) {
    printf("%14.6f  %12.6f  %12.6f\n", times[k], deltaRaCosDec[k], deltaDec[k]);
  }
}
