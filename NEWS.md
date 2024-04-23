# CptNonPar 0.2.0

# CptNonPar 0.2.0

* Added a new function multiscale.np.mojo, this function performs multiscale, multi-lag 
nonparametric change point detection using a set of bandwidths (as opposed 
to a single bandwidth as in np.mojo.mulilag).
* Updated parallelisation code in np.mojo to avoid unnecessary call to closeAllConnections(), thanks to 
Henrik Bengtsson for highlighting this.
* Changed uses of itemize in .Rd files to comply with CRAN check.

# CptNonPar 0.1.1

* Updated description field in Description file.
* Updated examples in np.mojo, np.mojo.multilag, and multilag.cpts.merge.

# CptNonPar 0.1.0

* Added a `NEWS.md` file to track changes to the package.
