# FiEstAS
Field Estimator for Arbitrary Spaces

The Field Estimator for Arbitrary Spaces (FiEstAS) computes the continuous probability density field underlying a given discrete data sample in multiple, non-commensurate dimensions.
The algorithm works by constructing a metric-independent tessellation of the data space based on a recursive binary splitting. Individual, data-driven bandwidths are assigned to each point, scaled so that a constant "mass" M<sub>0</sub> is enclosed.
Kernel density estimation may then be performed for different kernel shapes.
Further details may be found in:

- [Ascasibar & Binney (2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.356..872A/abstract), MNRAS 356, 872
- [Ascasibar (2010)](https://ui.adsabs.harvard.edu/abs/2010CoPhC.181.1438A/abstract), CoPhC, 181, 1438

*Please cite both articles if you use this software in your work.*

The code has not been properly tested or documented, but it should successfully compile with `make` under most systems and produce the executables `FiEstAS_ASCII`, `FiEstAS_Gadget`, and `randomData` in the `bin` directory. You may invoke them from the command line to see a short explanatory text describing usage. The easiest way to use the code is probably to save your data as an ASCII file and run

`FiEstAS_ASCII <your_file.txt>`

Please feel free to contact me for any help in using the code.
Any comments, suggestions, and/or (most importantly) bug reports are also welcome.

 Have fun :^)

  Yago

