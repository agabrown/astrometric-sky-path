# astrometric-sky-path
Use the SOFA library to calculate the coordinate directions of a source as a function of time, accounting for parallax and
proper motion.

The [`skypath.c`](./skypath.c) programme uses a few functions from the [SOFA](http://www.iausofa.org/) ANSI-C library to
calculate as a function of time the coordinate directions for a source with given astrometric catalogue parameters and radial
velocity. This results in the well-known helix- or wave-like motions of sources on the sky. The output of the programme is
directed to stdout, where it is picked up by the python script [`plotskypath.py`](./plotskypath.py) for making a plot of the
source path.

## Installation

The [SOFA ANSI-C library](http://www.iausofa.org/current_C.html#Downloads) should be installed. The makefile that is included
should be modified up to the point indicated in the file. Suggestions:

```INSTALL_DIR = $(HOME)```

```CFLAGF = -c -pedantic -Wall -W -O -fPIC```

In the latter line the `-fPIC` flag was added which is needed for the creation of a shared object library (`libsofa_c.so`).
If you do this remember to create the `.so` file with

```ld -G -z text -o libsofa_c.so *.o```

and copy it by hand to the installation folder (`$(HOME)/lib` in this example).

After installation of the SOFA libray do `make skypath` to create the `skypath` executable. Modify the [makefile](./makefile)
if needed.

## Usage

The `skypath` command can be invoked directly or through the python script. In both cases use the `-h` flag for usage
information. Note that the python script expects comma-separated lists of numbers between quotes as input for the
`--astrometry` or `--phaseSpace` command line arguments.
