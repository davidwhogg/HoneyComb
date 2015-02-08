# HoneyComb

Asteroseismology on the cheap.  Or really the expensive.

## Authors

* Ruth Angus (Oxford) (Harvard CfA)
* Dan Foreman-Mackey (NYU CCPP)
* David W. Hogg (NYU CCPP) (MPIA)

## License

Copyright 2014, 2015 the authors.
**HoneyComb** is open-source software released under the MIT License.
See the file `LICENSE` for details.

## Projects

* Stellar asteroseismology or parameter estimation without explicit line identification or peak identification in the periodogram.
* Replacement of the FFT or periodogram with a model that generates the data, including heterogenous, finite integration times on a heterogeneous grid, and heteroskedastic noise.
* Measurement of stellar oscillations at frequencies well above the classic "Nyquist" limit.

## Philosophy

In many data-analysis problems in the physical sciences,
the observations are (effectively) noisy projections of the world.
"The world" is itself a large set of continuous functions (fields)
of space, time, momentum, and so on.
The world has an essentially infinite number of degrees of freedom,
while the observations are always limited in number and signal-to-noise.
However, the things we *need* to know about the world are in many cases quite limited.
Good data analysis will infer what we need to infer
but avoid estimating all the details of all the continuous fields.
