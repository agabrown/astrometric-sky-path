/*
 * skypath_lib.c
 *
 * Author:        Anthony Brown
 * Last updated:  Oct 2018
 * First created: 23.08.2011
 *
 * Provide utility functions for the calculation of astrometric sky paths.
 */

#include "skypath_lib.h"

/*============================================================================
 *
 * usage()
 *
 * Input parameters: (from invocation of main(int argc, char **argv))
 * argv - command line arguments
 *
 *============================================================================
 */
void usage(char *argv[]) {
  fprintf(stderr,"Usage: %s ", argv[0]);
  fprintf(stderr,"-astrometry <comma separated list> -phaseSpace <comma separated list> -refepoch <reference epoch> -start <start epoch> -interval <time interval>\n");
  fprintf(stderr, "Use argument --help or -h for more help\n");
  exit(-1);
}

/*============================================================================
 *
 * usageLong()
 *
 * Input parameters: (from invocation of main(int argc, char **argv))
 * argv - command line arguments
 *
 *============================================================================
 */
void usageLong(char *argv[]) {
  fprintf(stderr,"Usage: %s [-args..]\n",argv[0]);
  fprintf(stderr,"  -astrometry   - Comma separated list of astrometric parameters\n");
  fprintf(stderr,"                  right ascension (degrees)\n");
  fprintf(stderr,"                  declination (degrees)\n");
  fprintf(stderr,"                  parallax (mas)\n");
  fprintf(stderr,"                  proper motion in RA  (mas/yr, mura*cos(delta))\n");
  fprintf(stderr,"                  proper motion in DEC (mas/yr)\n");
  fprintf(stderr,"                  radial Velocity (km/s)\n");
  fprintf(stderr,"  -phaseSpace   - Comma separated list of phase space parameters\n");
  fprintf(stderr,"                  X,Y,Z (ICRS, pc)\n");
  fprintf(stderr,"                  VX,VY,VZ (ICRS, km/s)\n");
  fprintf(stderr,"  -refepoch     - Reference epoch for astrometric parameters\n");
  fprintf(stderr,"  -start        - Start epoch for skypath calculation (Julian)\n");
  fprintf(stderr,"  -interval     - Time interval for skypath calculation (Julian years)\n");
  fprintf(stderr,"  -npoints      - Number of points to calculate between start and end epoch (default 1001)\n");
  fprintf(stderr,"  -h            - Show this help text\n");
  fprintf(stderr,"  --help        - Show this help text\n");
  exit(-1);
}

/*============================================================================
 *
 * parseArgs()
 *
 * Parse the command line arguments.
 *
 * Input parameters
 * argc - number of command line arguments (from invocation of main(int argc,
 *        char **argv))
 * argv - command line arguments
 * theStar - on return contains the star astrometric and phase space parameters.
 * obsEpochs - on return contains the observation time interval and reference epoch.
 *
 *============================================================================
 */
void parseArgs(int argc, char *argv[], int *npoints, star *theStar, obsEpochs *theObsEpochs) {
  int argCounter,k;
  int astrometryProvided=0;
  int phaseSpaceProvided=0;
  int astrometryOrPhaseSpace=1;
  char *token, *paramsString;
  double *starParams;
  star sofaStar;
  static const char DELIMITER[1]=",";
  static const int DEFAULT_NPOINTS=1001;

  *npoints = DEFAULT_NPOINTS;

  /* Check if help is asked for.*/
  if (argc>1) {
    for (argCounter=1;argCounter<argc;argCounter++) {
      if (strcmp(argv[argCounter], "--help") == 0 || strcmp(argv[argCounter], "-help") == 0
          || strcmp(argv[argCounter], "--h") == 0 || strcmp(argv[argCounter], "-h") == 0) {
        usageLong(argv);
      }
    }
  }

  initCalc(theStar, theObsEpochs);

  /* No command line arguments, run with default paramaters. */
  if (argc<2) return;

  starParams=dvector(6);

  /* Parse arguments */
  for (argCounter=1;argCounter<argc;argCounter++) {

    if (strcmp(argv[argCounter], "-astrometry") == 0) {
      astrometryProvided=1;
      astrometryOrPhaseSpace=1;
      if (argc < argCounter+2) usage(argv);
      paramsString = argv[++argCounter];
    }

    if (strcmp(argv[argCounter], "-phaseSpace") == 0) {
      phaseSpaceProvided=1;
      astrometryOrPhaseSpace=0;
      if (argc < argCounter+2) usage(argv);
      if (astrometryProvided==0) {
        paramsString = argv[++argCounter];
      }
    }

    if (strcmp(argv[argCounter], "-refepoch") == 0) {
      if (argc < argCounter+2) usage(argv);
      theObsEpochs->refEpJ=atof(argv[++argCounter]);
    }

    if (strcmp(argv[argCounter], "-start") == 0) {
      if (argc < argCounter+2) usage(argv);
      theObsEpochs->beginEpJ=atof(argv[++argCounter]);
    }

    if (strcmp(argv[argCounter], "-interval") == 0) {
      if (argc < argCounter+2) usage(argv);
      theObsEpochs->deltaEpJ=atof(argv[++argCounter]);
    }

    if (strcmp(argv[argCounter], "-npoints") == 0) {
      if (argc < argCounter+2) usage(argv);
      *npoints=atoi(argv[++argCounter]);
    }
  }

  if (astrometryProvided!=0 || phaseSpaceProvided!=0) {
      token = strtok(paramsString, DELIMITER);
      for (k=0; k<6; k++) {
          starParams[k]=atof(token);
          token = strtok(NULL, DELIMITER);
      }
  }

  if (astrometryProvided==1 && phaseSpaceProvided==1) {
    fprintf(stderr,"Both astrometry and phase coordinates provided; defaulting to astrometry input!\n");
    astrometryOrPhaseSpace=1;
  }

  if (astrometryOrPhaseSpace==1) {
      if (astrometryProvided==1) {
          sofaStar.rightAscension = starParams[0]/RADIAN_TO_DEGREE;
          sofaStar.declination = starParams[1]/RADIAN_TO_DEGREE;
          sofaStar.parallax = starParams[2]/1000.0;
          sofaStar.properMotionRa = starParams[3]/(1000.0*3600.0*RADIAN_TO_DEGREE)/cos(sofaStar.declination);
          sofaStar.properMotionDec = starParams[4]/(1000.0*3600.0*RADIAN_TO_DEGREE);
          sofaStar.radialVelocity = starParams[5];
      } else {
          myUnitsToSofa(theStar, &sofaStar);
      }
    iauStarpv(sofaStar.rightAscension, sofaStar.declination, sofaStar.properMotionRa,
        sofaStar.properMotionDec, sofaStar.parallax, sofaStar.radialVelocity,
        sofaStar.phaseSpaceCoordinates);
  } else {
    for (k=0; k<3; k++) {
      sofaStar.phaseSpaceCoordinates[0][k]=starParams[k]*PARSEC_TO_AU;
      sofaStar.phaseSpaceCoordinates[1][k]=starParams[k+3]*86400.0/AU_TO_KM;
    }
    iauPvstar(sofaStar.phaseSpaceCoordinates, &sofaStar.rightAscension, &sofaStar.declination,
        &sofaStar.properMotionRa, &sofaStar.properMotionDec, &sofaStar.parallax,
        &sofaStar.radialVelocity);
  }
  sofaUnitsToMine(&sofaStar, theStar);

  free_dvector(starParams);

}

/*===========================================================================
 *
 * initCalc()
 *
 * Initialize default star astrometric parameters, the corresponding phase space coordinates, and the
 * observation time interval.
 *
 * theStar - Variable of type "star"
 * theObsEpochs - Data structure holding the observation time interval and reference epoch.
 *
 *===========================================================================
 */
void initCalc(star *theStar, obsEpochs *theObsEpochs) {
  star sofaStar;
  sofaStar.rightAscension=30.0/RADIAN_TO_DEGREE;
  sofaStar.declination=40.0/RADIAN_TO_DEGREE;
  sofaStar.properMotionRa=40.0/(1000.0*3600.0*RADIAN_TO_DEGREE)/cos(sofaStar.declination);
  sofaStar.properMotionDec=30.0/(1000.0*3600.0*RADIAN_TO_DEGREE);
  sofaStar.parallax=100.0/1000.0;
  sofaStar.radialVelocity=30.0;
  iauStarpv(sofaStar.rightAscension, sofaStar.declination, sofaStar.properMotionRa,
      sofaStar.properMotionDec, sofaStar.parallax, sofaStar.radialVelocity,
      sofaStar.phaseSpaceCoordinates);
  sofaUnitsToMine(&sofaStar, theStar);
  theObsEpochs->refEpJ = 2017.0;
  theObsEpochs->beginEpJ = 2014.5;
  theObsEpochs->deltaEpJ = 5.0;
}

/*===========================================================================
 *
 * myUnitsToSofa()
 *
 * Convert the struct star's parameters from the units used by default in the struct to the units used by
 * SOFA.
 *
 * myStar - struct star with default units (input)
 * sofaStar - struct star with SOFA units (output)
 *
 *===========================================================================
 */
void myUnitsToSofa(star *myStar, star *sofaStar) {
  int k;
  sofaStar->rightAscension = myStar->rightAscension;
  sofaStar->declination = myStar->declination;
  sofaStar->properMotionRa =
    myStar->properMotionRa/(1000.0*3600.0*RADIAN_TO_DEGREE)/cos(myStar->declination);
  sofaStar->properMotionDec = myStar->properMotionDec/(1000.0*3600.0*RADIAN_TO_DEGREE);
  sofaStar->parallax = myStar->parallax/1000.0;
  sofaStar->radialVelocity = myStar->radialVelocity;
  for (k=0; k<3; k++) {
    sofaStar->phaseSpaceCoordinates[0][k]=myStar->phaseSpaceCoordinates[0][k]*PARSEC_TO_AU;
    sofaStar->phaseSpaceCoordinates[1][k]=myStar->phaseSpaceCoordinates[1][k]*86400.0/AU_TO_KM;
  }
}

/*===========================================================================
 *
 * sofaUnitsToMine()
 *
 * Convert the struct star's parameters from the units used by default in the struct to the units used by
 * SOFA.
 *
 * sofaStar - struct star with SOFA units (input)
 * myStar - struct star with default units (output)
 *
 *===========================================================================
 */
void sofaUnitsToMine(star *sofaStar, star *myStar) {
  int k;
  myStar->rightAscension = sofaStar->rightAscension;
  myStar->declination = sofaStar->declination;
  myStar->properMotionRa =
    sofaStar->properMotionRa*(1000.0*3600.0*RADIAN_TO_DEGREE)*cos(sofaStar->declination);
  myStar->properMotionDec = sofaStar->properMotionDec*(1000.0*3600.0*RADIAN_TO_DEGREE);
  myStar->parallax = sofaStar->parallax*1000.0;
  myStar->radialVelocity = sofaStar->radialVelocity;
  for (k=0; k<3; k++) {
    myStar->phaseSpaceCoordinates[0][k]=sofaStar->phaseSpaceCoordinates[0][k]/PARSEC_TO_AU;
    myStar->phaseSpaceCoordinates[1][k]=sofaStar->phaseSpaceCoordinates[1][k]/86400.0*AU_TO_KM;
  }
}

/*===========================================================================
 *
 * calcSkyPath()
 *
 * Calculate the astrometric path on the sky.
 *
 * theObsEpochs - Observation time interval and reference epoch
 * npoints - number of points on path.
 * sofaStar - struct star IN SOFA UNITS
 * times - output time array (Julian Epoch)
 * deltaRaCosDec - output array with values of delta(RA)*cos(DEC)
 * deltaDec - output array with values of delta(DEC)
 *
 *===========================================================================
 */
void calcSkyPath(obsEpochs *theObsEpochs, int npoints, star *sofaStar, double *times, double*
    deltaRaCosDec, double *deltaDec){

  int k;
  double heliocentricEarth[2][3], barycentricEarth[2][3], pvMatrix[2][3];
  double ra, dec;
  double raZero, decZero, pmraZero, pmdecZero, parallaxZero, rvZero;
  double refModifiedJD;
  double zeropointJD, modifiedJD;
  double coordinateDirection[3];

  for (k=0; k<3; k++) {
    pvMatrix[0][k]=sofaStar->phaseSpaceCoordinates[0][k];
    pvMatrix[1][k]=sofaStar->phaseSpaceCoordinates[1][k];
  }

  iauPvstar(pvMatrix, &raZero, &decZero, &pmraZero, &pmdecZero, &parallaxZero, &rvZero);

  iauEpj2jd(theObsEpochs->refEpJ, &zeropointJD, &refModifiedJD);

  for (k=0; k< npoints; k++) {
    times[k] = theObsEpochs->beginEpJ+k*(theObsEpochs->deltaEpJ/(npoints-1));
    /*
     * Calculate the earth's position in ICRS coordinates centred on the Solar
     * system barycentre.
     */
    iauEpj2jd(times[k], &zeropointJD, &modifiedJD);
    iauEpv00(zeropointJD, modifiedJD, heliocentricEarth, barycentricEarth);
    /*
     * Calculate the coordinate direction to the star accounting only for the proper motion and parallax.
     * This gives the path on the sky as (alpha(t), delta(t)) without accounting for all the other
     * effects that determine the actually observed direction at epoch t.
     */
    iauPmpx(sofaStar->rightAscension, sofaStar->declination, sofaStar->properMotionRa,
            sofaStar->properMotionDec, sofaStar->parallax, sofaStar->radialVelocity,
            times[k]-theObsEpochs->refEpJ, barycentricEarth[0], coordinateDirection);
    iauC2s(coordinateDirection, &ra, &dec);
    /*
     * Calculate the sky path as (delta_ra*cos(dec), delta_dec) in mas and store in an output array. Take
     * care of the cases where the star crosses the RA=0 or RA=180 degrees line.
     */
    deltaRaCosDec[k]=ra-raZero;
    if (deltaRaCosDec[k] > PI) deltaRaCosDec[k]=deltaRaCosDec[k]-2*PI;
    if (deltaRaCosDec[k] < -PI) deltaRaCosDec[k]=deltaRaCosDec[k]+2*PI;
    deltaDec[k]=dec-decZero;
    deltaRaCosDec[k]=deltaRaCosDec[k]*RADIAN_TO_DEGREE*3600.0*1000.0*cos(dec);
    deltaDec[k]=deltaDec[k]*RADIAN_TO_DEGREE*3600.0*1000.0;
  }
}
