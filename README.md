# astrometric-sky-path
Use the SOFA library to calculate the coordinate directions of a source as a function of time, accounting for parallax and
proper motion.

_Anthony G.A. Brown <brown@strw.leidenuniv.nl>_

The [`skypath.c`](./skypath.c) programme uses a few functions from the [SOFA](http://www.iausofa.org/) ANSI-C library to
calculate as a function of time the coordinate directions, as seen from the observers position, for a source with given
astrometric catalogue parameters and radial velocity. This results in the well-known helix- or wave-like motions of sources
on the sky. The output of the programme is directed to stdout, where it is picked up by the python script
[`plotskypath.py`](./plotskypath.py) for making a plot of the source path. The plot of the source path is
done using so-called local plane coordinates which are defined by equation 1.2.22 in section 1.2 of the
[Hipparcos documentation](https://ui.adsabs.harvard.edu/#abs/1997ESASP1200.....E/abstract) (see [this
link](https://www.cosmos.esa.int/documents/532822/552851/vol1_all.pdf/99adf6e3-6893-4824-8fc2-8d3c9cbba2b5)
for the PDF document). These local plane coordinates are also calculated by the `skypath.c` programme.

This primitive Python to C interfacing was a lazy option to quickly make a plot. The main purpose of this exercise was to
learn more about using SOFA.

__NOTE__ that the SOFA library is also used in [Astropy](http://www.astropy.org/) (in derived form, see [Astropy
ERFA](https://github.com/astropy/astropy/tree/master/cextern/erfa)), but the calculation of coordinate directions including
parallax and proper motion effects is not yet implemented as far as I can tell.

## Installation

The [SOFA ANSI-C library](http://www.iausofa.org/current_C.html#Downloads) should be installed. The makefile that is included
should be modified up to the point indicated in the file. Suggestions:

```INSTALL_DIR = $(HOME)```

```CFLAGF = -c -pedantic -Wall -W -O -fPIC```

In the latter line the `-fPIC` flag was added which is needed for the creation of a shared object library (`libsofa_c.so`).
If you do this remember to create the `.so` file with

```ld -G -z text -o libsofa_c.so *.o```

and copy it by hand to the installation folder (`$(HOME)/lib` in this example).

After installation of the SOFA library do `make skypath` to create the `skypath` executable. Modify the [makefile](./makefile)
if needed.

## Usage

The `skypath` command can be invoked directly or through the python script. In both cases use the `-h` flag for usage
information. Note that the python script expects comma-separated lists of numbers between quotes as input for the
`--astrometry` or `--phaseSpace` command line arguments.
