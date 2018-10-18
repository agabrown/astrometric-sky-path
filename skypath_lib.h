/*
 * Header file for skypath_lib.c
 *
 * Anthony Brown 23 Aug 2011 - Oct 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sofa.h"
#include "sofam.h"
#include "vecmat_utils.h"
#include "astrometry.h"

// Rename a number of the SOFA constants for clarity
static const double RADIAN_TO_DEGREE = DR2D;
static const double PARSEC_TO_AU = DR2AS;
static const double AU_TO_KM = DAU/1000.0;
static const double PI = DPI;

/* Function prototypes */

void usage(char *argv[]);

void usageLong(char *argv[]);

void parseArgs(int argc, char *argv[], int *npoints, star *theStar, obsEpochs *theObsEpochs);

void initCalc(star *theStar, obsEpochs *theObsEpochs);

void myUnitsToSofa(star *myStar, star *sofaStar);

void sofaUnitsToMine(star *sofaStar, star *myStar);

void calcSkyPath(obsEpochs *theObsEpochs, int npoints, star *sofaStar, double *times, double* deltaRaCosDec, double
        *deltaDec);
