# cppNumDeriv
Methods for calculating numerical first and second order derivatives. This header-only lib is a C++ 11 (partial) porting of its R counterpart [numDeriv-package](https://cran.r-project.org/web/packages/numDeriv) by Paul Gilbert and Ravi Varadhan.
Accurate calculations are done using 'Richardson''s' extrapolation. A simple difference method is also provided. Simple difference is (usually) less accurate but is much quicker.
## The main functions are
* __gradient__   to calculate the gradient (first derivative) of a scalar real valued function with real n-vector argument. If method is "simple", the calculation is done using a simple epsilon difference. If method is "Richardson", the calculation is done by Richardson’s extrapolation (see e.g.  Linfield and Penny, 1989, or Fornberg and Sloan, 1994.) 
* __hessian__   to calculate the Hessian (second derivative) of a scalar real valued function with real n-vector argument. The "Richardson" method uses genD to extract the second derivative.
* __genD__   to calculate the gradient and second derivative of a real m-vector valued function with real n-vector argument. genD generates a matrix of function derivative information. The derivatives are calculated numerically using Richardson improvement. The "Richardson" method calculates a numerical approximation of the first and second derivatives of func at the point x.

## References
See [numDeriv-package](https://cran.r-project.org/web/packages/numDeriv/numDeriv.pdf) for a detailed overview.

## References
cppNumDeriv depends on the [Eigen](http://eigen.tuxfamily.org/index.php) template library. Visit this link for a detailed overview.
