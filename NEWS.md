# CptNonPar 0.3.0

# Changes in version 0.3.0

* The associated paper is now published in Biometrika: see 
<doi:10.1093/biomet/asaf024> for full details.
* Updated `multiscale.np.mojo()` function so that returned cpts are given in time order.
* The p-values returned by `np.mojo()`, `np.mojo.multilag()`, and `mulsticale.np.mojo()` 
have been replaced by importance scores: when the bootstrap is used, these are essentially one minus the p-value. 
* Fixed error when using the manual threshold and the `np.mojo.multilag()` function, thanks to 
Chuanyang Zhang for spotting this.
* You can now use different manual thresholds for different lags: for use 
in `np.mojo.multilag()`, you can supply a vector, whilst for `mulsticale.np.mojo()`, you can supply a list of vectors.
* Added package level documentation: see `?CptNonPar`.

# Changes in version 0.2.1

* Updated link to the paper in description to comply with CRAN check.

# Changes in version 0.2.0

* Added a new function `multiscale.np.mojo()`, this function performs multiscale, multi-lag 
nonparametric change point detection using a set of bandwidths (as opposed 
to a single bandwidth as in `np.mojo.mulilag()`).
* Updated parallelisation code in `np.mojo()` to avoid unnecessary call to `closeAllConnections()`, thanks to 
Henrik Bengtsson for highlighting this.
* Changed uses of itemize in .Rd files to comply with CRAN check.

# Changes in version 0.1.1

* Updated description field in Description file.
* Updated examples in np.mojo, np.mojo.multilag, and multilag.cpts.merge.

# Changes in version 0.1.0

* Added a `NEWS.md` file to track changes to the package.
