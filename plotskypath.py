"""
Plot the astrometric path of stars on the sky. That is plot the coordinate direction including the proper
motion and parallax effect as a function of time.

Use the C-program skypath.c to draw astrometric paths of stars on the sky. Requires the SOFA library:
http://www.iausofa.org/ (ANSI-C version).

Anthony Brown 2011-2018
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse

# Configure matplotlib
rc('text', usetex=True)
rc('font', family='serif', size=18)
rc('xtick.major', size='10')
rc('xtick.minor', size='5')
#rc('xtick', direction='out')
rc('ytick.major', size='10')
rc('ytick.minor', size='5')
#rc('ytick', direction='out')
rc('lines', linewidth=2)
rc('axes', linewidth=2)


def make_plot(args):
    """
    Calculate the astrometric path on the sky and plot the result.

    Parameters
    ----------

    args - Command line arguments.

    Returns
    -------

    Nothing
    """

    skyPathCommand=['./skypath']
    refEpoch = np.float(args['refEpochString'])
    startEpoch = np.float(args['startEpochString'])
    npoints = args['npoints']

    if (args['parEllipse']):
        skyPathCommandParallaxOnly=['./skypath']
    if (args['astrometryString']!=None):
        skyPathCommand.append('-astrometry')
        skyPathCommand.append(args['astrometryString'])
    if (args['phaseSpaceString']!=None):
        skyPathCommand.append('-phaseSpace')
        skyPathCommand.append(args['phaseSpaceString'])
    if (args['refEpochString']!=None):
        skyPathCommand.append('-refepoch')
        skyPathCommand.append('{0}'.format(args['refEpochString']))
    if (args['startEpochString']!=None):
        skyPathCommand.append('-start')
        skyPathCommand.append('{0}'.format(args['startEpochString']))
    if (args['timeIntervalString']!=None):
        skyPathCommand.append('-interval')
        skyPathCommand.append('{0}'.format(args['timeIntervalString']))
    if (args['npoints']!=None):
        skyPathCommand.append('-npoints')
        skyPathCommand.append('{0}'.format(args['npoints']))

    if (args['parEllipse']):
        # TODO Better defaults handling in this case.
        if (args['astrometryString']==None and args['phaseSpaceString']==None):
            args['astrometryString'] = "30,40,100,40,30,30"
        if (args['astrometryString']!=None):
            skyPathCommandParallaxOnly=['./skypath']
            skyPathCommandParallaxOnly.append('-astrometry')
            alpha,delta,parlx,mualpha,mudelta,vrad=args['astrometryString'].split(",")
            parOnlyString='{0},{1},{2},0,0,0'.format(alpha,delta,parlx)
            skyPathCommandParallaxOnly.append(parOnlyString)
            skyPathCommandParallaxOnly.append('-interval')
            skyPathCommandParallaxOnly.append('1.0')
        if (args['refEpochString']!=None):
            skyPathCommandParallaxOnly.append('-refepoch')
            skyPathCommandParallaxOnly.append('{0}'.format(args['refEpochString']))
        if (args['phaseSpaceString']!=None):
            skyPathCommandParallaxOnly.append('-astrometry')
            x,y,z,vx,vy,vz=args['astrometryString'].split(",")
            parOnlyString='{0},{1},{2}.0.0,0.0,0.0'.format(x,y,z)
            skyPathCommandParallaxOnly.append(parOnlyString)
            skyPathCommandParallaxOnly.append('-interval')
            skyPathCommandParallaxOnly.append('1.0')
        if (args['startEpochString']!=None):
            skyPathCommandParallaxOnly.append('-start')
            skyPathCommandParallaxOnly.append('{0}'.format(args['startEpochString']))

    result=subprocess.run(skyPathCommand, stdout=subprocess.PIPE)
    skyPath=result.stdout.splitlines()
    times=np.empty(len(skyPath))
    deltaRaCosDec=np.empty(len(skyPath))
    deltaDec=np.empty(len(skyPath))
    for i in range(len(skyPath)):
        times[i], deltaRaCosDec[i], deltaDec[i] = skyPath[i].split()

    if (args['parEllipse']):
        resultB=subprocess.run(skyPathCommandParallaxOnly, stdout=subprocess.PIPE)
        skyPathB=resultB.stdout.splitlines()
        timesB=np.empty(len(skyPathB))
        deltaRaCosDecParOnly=np.empty(len(skyPathB))
        deltaDecParOnly=np.empty(len(skyPathB))
        for i in range(len(skyPathB)):
            timesB[i], deltaRaCosDecParOnly[i], deltaDecParOnly[i] = skyPathB[i].split()

    fig=plt.figure(figsize=(8,8))
    if (args['parEllipse']):
        plt.plot(deltaRaCosDecParOnly, deltaDecParOnly,'k--',alpha=0.5)
    if args['plotDots']:
        plt.plot(deltaRaCosDec, deltaDec, 'o')
    else:
        plt.plot(deltaRaCosDec, deltaDec)
    plt.scatter(deltaRaCosDec[0], deltaDec[0], c='r', marker='o', s=50)
    indref = np.searchsorted(times, refEpoch, side='right')-1
    plt.scatter(deltaRaCosDec[indref], deltaDec[indref], c='r', marker='+', s=50)
    plt.scatter(deltaRaCosDec[-1], deltaDec[-1], c='r', marker='^', s=50)

    plt.xlabel("$\\Delta\\alpha\\cos\\delta$ [mas]")
    plt.ylabel("$\\Delta\\delta$ [mas]")
    plt.grid()
    if (args['axisLimits']!=None):
        plt.xlim(args['axisLimits'])
        plt.ylim(args['axisLimits'])

    basename = 'skyPathIcrs'
    if args['pdfOutput']:
        plt.savefig(basename+'.pdf')
    elif args['pngOutput']:
        plt.savefig(basename+'.png')
    else:
        plt.show()

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser("Draw astrometric paths of stars on the sky.")
    parser.add_argument("--astrometry", dest="astrometryString",
    help="""Comma-separated list of astrometric parameters
          (alpha, delta, parallax, mu_alpha*cos(delta), mu_delta, Vrad)
          [deg, deg, mas, mas/yr, mas/yr, km/s]""")
    parser.add_argument("--phaseSpace", dest="phaseSpaceString",
            help="""Comma-separated list of phase space coordinates (X, Y, Z, Vx, Vy, Vz) [pc, pc, pc,
            km/s, km/s, km/s]""")
    parser.add_argument("--ref", dest="refEpochString", help="Reference epoch (Julian years)",
            default=2017.0)
    parser.add_argument("--start", dest="startEpochString", help="Start epoch (Julian years)",
            default=2014.5)
    parser.add_argument("--interval", dest="timeIntervalString", help="time interval (Julian years)",
            default=5.0)
    parser.add_argument("--plotLimits", dest="axisLimits", type=float, nargs=2, help="list of plot axis limits (low high) for both x and y")
    parser.add_argument("--parEllipse", action="store_true", dest="parEllipse", help="Show parallax ellipse")
    parser.add_argument("--npoints", help="""Number of points to calculate between start and end
            epoch""", type=int, default=1001)
    parser.add_argument("-d", action="store_true", dest="plotDots", help="""Plot individual points instead
            of a continuous line.""")
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
    parser.add_argument("-c", action="store_true", dest="colourFigure", help="Make colour plot")
    parser.add_argument("-t", action="store_true", dest="forTalk",  help="make version for presentations")
    args = vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    make_plot(args)
