# Interstellar gas analysis

Analysis work from graduate Gaseous Astrophysics (TCU, Fall 2022). The through-line is
measuring the properties of interstellar gas from spectra — how much gas there is, how fast it
is moving, and what is ionizing it — across both emission and absorption, and in both single
spectra and full integral-field datacubes.

Everything here is Python: NumPy, SciPy, AstroPy, Matplotlib.

## The analyses

**HI 21cm Gaussian decomposition** — `Pre-Projects/Pre-Project 1/Code/Gaussian_Code.py`

Decomposes neutral-hydrogen emission spectra from the HI4PI survey into overlapping Gaussian
velocity components along a strip of the sky, then converts each component to a hydrogen column
density. The interesting part is that the number of components is not the same from sightline to
sightline — blended emission has to be separated into physically distinct clouds before the
column densities mean anything, and the component count and initial parameter guesses are set
per sightline. Column densities come from integrating each fitted Gaussian analytically with the
error function over a finite velocity range rather than summing pixels, so the result does not
depend on where the integration window is clipped. Outputs are the velocity centroid, line
width, and column density as functions of Galactic longitude.

```bash
python Gaussian_Code.py --data-dir path/to/HI4PI --savefigs
```

**Emission-line fitting and AGN classification in IFU datacubes** —
`Pre-Projects/Pre_Project_2/` and `Projects/Project_2/`

Works with MaNGA integral-field spectroscopy, where every spatial pixel carries its own
spectrum. Fits multi-component Gaussians to the nebular emission lines, then uses the
line-ratio diagnostics (BPT diagrams in both the [N II] and [S II] flavors) to classify each
region by what is exciting the gas — ordinary star formation, or an active galactic nucleus.
The classification boundaries are implemented directly rather than eyeballed off a plot.

**UV absorption-line column densities** — `Pre-Projects/Pre-Project_3/` and `Projects/Project_3/`

The absorption side of the same question, on HST/COS far-UV spectra. Converts heliocentric to
Local Standard of Rest velocities, fits and divides out the stellar continuum with polynomials,
integrates apparent optical depth with Simpson's rule to get column densities, and identifies
transitions against a line list. `Projects/Project_3` adds a **Monte Carlo** treatment to
propagate measurement uncertainty into the derived columns instead of relying on a single
best-fit value. Outputs include the velocity-aligned column-density and flux plotstacks
(`Col_Dens_Plotstack.svg`, `Col_Flux_Plotstack.svg`).

**HST observing proposal** — `HST_Proposal/`

A full proposal for UV absorption spectroscopy of M31's Fermi Bubbles, which is target selection
as a data problem: query GALEX for UV-bright background sources (ten catalog searches, in
`GALEX Search Files/`), cross-check which are actually bright enough to yield usable S/N in
reasonable exposure time, and pull LAB survey HI spectra along the surviving sightlines to
confirm there is foreground gas to measure. `Project 1 Jupyter Notebook.ipynb` also carries the
gas mass estimation, and the Galactic All-Sky Survey HI cube is used here to check the
foreground. `Projects/Project_1/Code/HI_Fits_Analysis.py` is a short FITS-loading stub for that
survey rather than an analysis in its own right.

## Running any of this

The scripts take their data directories as arguments or environment variables and default to a
`Data/` directory beside the code. Survey data is not committed — HI4PI, GASS, and LAB HI cubes
and MaNGA datacubes are large and available from their own archives, and HST/COS spectra come
from [MAST](https://mast.stsci.edu). The committed figures and SVG plotstacks show the outputs.

`Gaussian_Code.py` additionally uses [gaussdecomp](https://github.com/dnidever/gaussdecomp).
