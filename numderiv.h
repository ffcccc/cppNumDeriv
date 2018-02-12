#ifndef ___numDerviHeader___
#define ___numDerviHeader___

#include <iostream>
#include <tuple>
#include <algorithm>
#include <functional>
#include "Eigen/Core"
#include "EigenIter/eigenBeginEnd.h"

//using namespace Eigen;

typedef Eigen::ArrayXd ValArr;
typedef Eigen::ArrayXXd ValMat;

// helper: print error string & return errorcode
int stop(const std::string& errmessage, int errcode = 0)
{
    std::cerr << errmessage;
    return errcode;
}

// default NA value
#define NA 1000000.0

//	default method.args in numderiv R package:
//	(eps = 1e-4, d = 0.0001, zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2, show.details = FALSE)

// method.args labels:
#define _EPS 0
#define _D 1
#define _ZERO_TOL 2
#define _R 3
#define _V 4
#define _SHOW_DETAILS 5

// method.args tuple
typedef std::tuple<double, double, double, int, double, bool> TArgs;
// TArgs args = std::make_tuple(0.0001, 0.0001, sqrt(.Machine$double.eps / 7e-7), 4, 2, false);
// TArgs args = std::make_tuple(0.0001, 0.0001, sqrt(std::numeric_limits<T>::epsilon() / 7e-7), 4, 2, false);
// TArgs args = std::make_tuple(0.0001, 0.0001, 1.781e-5, 4, 2, false);


// easy class to encapsulate F: Rnx --> Rny
// assert DimY = 1
//
class TFunc
{
public:
    // take the field to sort by in the constructor
    TFunc(int nx, int ny, std::function<double(const ValArr&)> func_par) : _dimX(nx) , _dimY(ny)
    {
		myfunc = func_par;
    }

    // case 1,2
    virtual void operator()(const ValArr& x, double& f) {
		//f = x.sin().sum();
		f = myfunc(x);
    };

    inline int DimX()
    {
	return _dimX;
    };

    inline int DimY()
    {
	return _dimY;
    };
    void SetDimX(int _dim)
    {
	this->_dimX = _dim;
    }
	void SetDimY(int _dim)
    {
	this->_dimY = _dim;
    }

private:
	// functor
    std::function<double(const ValArr&)> myfunc;
	// domain, codomain dims
    int _dimX;
    int _dimY;
};

//############################################################################
//
//#    functions for gradient calculation case 2
//
//############################################################################

//  # case 1/ scalar arg, scalar result (case 2/ or 3/ code should work)
//  # case 2/ vector arg, scalar result (same as special case jacobian)
//  # case 3/ vector arg, vector result (of same length, really 1/ applied multiple times))

//   #   gradient case 1 and 2 are special cases of jacobian:
//   #   with a scalar rather than vector valued function.
//   #   Case 3 differs only because of the interpretation
//   #   that the vector result is a scalar function applied to each argument, and the
//   #   thus the result has the same length as the argument.
//   #   Case 3 is NOT implemented in this version. Only scalar results !
//   #   Jacobian is NOT implemented in this version.
//   #   i) The code of gradient could be consolidated to be a special case of jacobian.
//   #   ii)The code of gradient could be used to implement a Jacobian function
//
//    C++ porting by Fabio Rosa			
//    modified by Paul Gilbert
//    from original code by Xingqiao Liu

int gradient(TFunc func,
    const ValArr& x,
    ValArr& res,
    const std::string& method,
    const ValArr& side_par,
    const TArgs& method_args) // ,...
{
	int ny = func.DimY();
	int n = func.DimX();
    //int n = x.size(); // number of variables in argument
    double	epsilon	= std::get<_EPS>(method_args);
    int 	r		= std::get<_R>(method_args);

    ValArr side(side_par);
    // std::valarray<double> x(x_par.data(), n);

    if(side.size() == 0) {
		side.resize(n);
		side.setConstant(n, NA);
    } else {
		if(n != side.size())
			return stop("Non-NULL argument 'side' should have the same length as x");

		if(std::any_of(begin(side), end(side), [](double i) { return (i == NA) || (fabs(i) != 1.0); }))
		return stop("Non-NULL argument 'side' should have values NA, +1, or -1.");
    }

    //    bool case1  = (n == ny == 1);
    //    bool case2  = (n > 1)   && (ny == 1);
    //    bool case3  = (n == ny) && (n > 1);
    bool case1or3 = (n == ny);

    if((1 != ny) && !case1or3)
		return stop("grad assumes a scalar valued function.");

	res.resize(n);

	// #very simple numerical approximation
    if(method == "simple") {
		
		// side[is.na(side)] <- 1
		std::replace_if(begin(side), end(side), [](double i) { return (i == NA); }, 1.0);

		// eps<-rep(args$eps, n) * side
		//ValArr eps = side * epsilon;

		// #now case 2
		double f1,f2;
		func(x, f1);

		// df <-rep(NA, n)  we use a[0,] instead
		ValArr dx = x;
		for(int i = 0; i < n; i++) { // in 1 : n)
			dx[i] = x[i] + side[i] * epsilon;	//eps[i];
			func(dx, f2);
			dx[i] = x[i];
			
			res[i] = (f2 - f1) / epsilon;
		}

		// return ( (n1 - func(x, ...)) / eps)
		return 1;
    }
    else if(method == "Richardson") {

		double tol = std::get<_ZERO_TOL>(method_args);
		double d = std::get<_D>(method_args);
		double v = std::get<_V>(method_args);
		bool show_details = std::get<_SHOW_DETAILS>(method_args);

		// #first order derivatives are stored in the matrix a[k, i],
		// #where the indexing variables k for rows(1 to r), i for columns(1 to n),
		// #r is the number of iterations, and n is the number of variables.
		
		std::vector<ValArr> y(r);
		for(int k = 0; k < r; k++) {
			y[k].resize(n);
			y[k].setZero(); // apply(clear);
		}

		// h = abs(x * d);
		ValArr h = abs(x * d);
		//h *= d;
		for(int i = 0; i < n; i++) {
			if(fabs(x[i]) < tol) {
				h[i] += epsilon;
			}
		}
		ValArr ph(n);
		ValArr mh(n);
		ValArr dx(x);
		double f1, f2, tmp;

		// #successively reduce h
		ph = h;
		mh = h;
		for(int i = 0; i < n; i++) {
			if(side[i] == 1.0) {
				ph[i] = 2.0 * h[i];
				mh[i] = 0.0;
			} else if(side[i] == -1.0) {
				mh[i] = 2.0 * h[i];
				ph[i] = 0.0;
			}
		}
		for(int k = 0; k < r; k++) {
			for(int i = 0; i < n; i++) {
				
				if((k != 0) && (fabs(y[k - 1][i]) < 0.00000000000000000001)) {
					// #some func are unstable near zero
					y[k][i] = 0.0; 
				} else {
					dx[i] = x[i] + ph[i];
					func(dx, f1);

					dx[i] = x[i] - mh[i];
					func(dx, f2);

					dx[i] = x[i];
					tmp = (f1 - f2) / (2.0 * h[i]);
					y[k][i] = tmp;
				}
			}

			h /= v; // Reduced h by 1/v.
			ph /= v;
			mh /= v;
		}

		if(show_details) {
		}
		// #-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
		// # 1 Applying Richardson Extrapolation to improve the accuracy of
		// #the first and second order derivatives.The algorithm as follows:
		// #
		// #-- For each column of the derivative matrix a,
		// #say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
		// #new sequence of approximations B1, B2, ..., Br used the formula
		// #
		// #B(i) = (A(i + 1) * 4 ^ m - A(i)) / (4 ^ m - 1), i = 1, 2, ..., r - m
		// #
		// #N.B.This formula assumes v = 2.
		// #
		// #-- Initially m is taken as 1 and then the process is repeated
		// #restarting with the latest improved values and increasing the
		// #value of m by one each until m equals r - 1
		// #
		// # 2 Display the improved derivatives for each
		// #m from 1 to r - 1 if the argument show.details = T.
		// #
		// # 3 Return the final improved derivative vector.
		// #-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
		// -
		int rmax = r - 2;
		for(int m = 1; m <= r - 1; m++) {
			for(int i = rmax; i >= 0; i--) {
				double p4m = pow(4.0, m);
				y[i] = (y[i + 1] * p4m - y[i]) / (p4m - 1.0);
				// a = (a [2:(r+1-m), , drop = FALSE] * (4 ^ m) - a [1:(r - m), , drop = FALSE]) / (4 ^ m - 1)

				if(show_details) {
				}
			}
			rmax--;
		}

		res = y[0];
		return 1;
    } 
	else
		return stop("indicated method not supported.");
}

ValArr gradient(TFunc func,
    const ValArr& x,
    const std::string& method,
    const ValArr& side_par,
    const TArgs& method_args)
{
	int n = func.DimX();
	ValArr res(n);
	gradient(func, x, res, method, side_par, method_args);
	return res;
}

ValArr genD(TFunc func, const ValArr& x, const std::string& method, const TArgs& method_args)
{
	//    C++ porting by Fabio Rosa				(Jan, 2018)
    //    additional cleanup by Paul Gilbert 	(March, 2006)
    //    modified substantially by Paul Gilbert(May, 1992)
    //    from original code by Xingqiao Liu.	(May, 1991).
    //  
    //    The R original function was not optimized for speed, but organized in
    //    the same way it could be (was) implemented in C, to facilitate checking.
    //    The current C++ 11 version adopts the optimized Eigen library, but its architecture
	//    still needs more attention to speed
	
	//  
    //    v  reduction factor for Richardson iterations. This could
    //    be a parameter but the way the formula is coded it is assumed to be 2.

    if(method != "Richardson")
		stop("method not implemented.");

	double	d = std::get<_D>(method_args);
	double	v = std::get<_V>(method_args);
	int		r = std::get<_R>(method_args);
	int n = func.DimX();
	double epsilon = std::get<_EPS>(method_args);
	double tol = std::get<_ZERO_TOL>(method_args);

	if(v != 2)
		stop("The current code assumes v is 2 (the default).");

    // #f0 is the value of the function at x.
    double f0(0.0);
    func(x, f0);

    // n <- length(x)  #  number of parameters (theta)
    // h0 <- abs(d*x) + args$eps * (abs(x) < args$zero.tol)
    ValArr h0 = abs(x * d);
    for(int i = 0; i < n; i++) {
		if(fabs(x[i]) < tol) {
			h0[i] += epsilon;
		}
    }

    ValArr res((n * (n + 3)) / 2);
	res.setZero();
	
    // length(f0) is the dim of the sample space
    // (n*(n + 3))/2 is the number of columns of matrix D.( first
    //	der. & lower triangle of Hessian)
    ValArr Daprox(r);	Daprox.setZero();
    ValArr Hdiag(n);	Hdiag.setZero();
    ValArr Haprox(r);	Haprox.setZero();
    ValArr h(n);
    ValArr dx(n);
	
    // each parameter  - first deriv. & hessian diagonal
    for(int i = 0; i < n; i++) {
		h = h0;
		dx = x;
		double f1, f2;
		for(int k = 0; k < r; k++) { // successively reduce h
			dx[i] = x[i] + h[i];
			func(dx, f1);
			dx[i] = x[i] - h[i];
			func(dx, f2);
			dx[i] = x[i];
			// double f1; func(x+(i==(1:n))*h, f1)
			// double f2; func(x-(i==(1:n))*h, f2)

			Daprox[k] = (f1 - f2) / (2 * h[i]);            // F'(i)
			Haprox[k] = (f1 - 2 * f0 + f2) / pow(h[i], 2); // F''(i,i) hessian diagonal
			h /= v;                                        // Reduced h by 1/v.
		}

		int rmax = r - 2;
		for(int m = 1; m <= r - 1; m++) {
			for(int k = rmax; k >= 0; k--) {
				double p4m = pow(4, m);
				Daprox[k] = (Daprox[k + 1] * p4m - Daprox[k]) / (p4m - 1);
				Haprox[k] = (Haprox[k + 1] * p4m - Haprox[k]) / (p4m - 1);
			}
			rmax--;
		}
		res[i] = Daprox[0];
		Hdiag[i] = Haprox[0];
    }

	h.setZero();
	dx.setZero();
	Daprox.setZero();
    int u = n;

    for(int i=0; i<n; i++) { // 2nd derivative  - do lower half of hessian only
		
		for(int j=0; j<=i; j++) {

			if(i == j) {
				res[u] = Hdiag[i];
			}
			else {
				h = h0;
				dx = x;
				double f1, f2;
				for(int k = 0; k < r; k++) { // successively reduce h
					dx[i] = x[i] + h[i];
					dx[j] = x[j] + h[j];
					func(dx, f1);
					dx[i] = x[i] - h[i];
					dx[j] = x[j] - h[j];
					func(dx, f2);
					dx[i] = x[i];
					dx[j] = x[j];
					Daprox[k] = (f1 - 2 * f0 + f2 - Hdiag[i] * pow(h[i], 2) - Hdiag[j] * pow(h[j], 2)) / (2 * h[i] * h[j]); // F''(i,j)
					h /= v;                // Reduced h by 1/v.
				}

				int rmax = r - 2;
				for(int m = 1; m <= r - 1; m++) {
					for(int k = rmax; k >= 0; k--) {
						double p4m = pow(4, m);
						Daprox[k] = (Daprox[k + 1] * p4m - Daprox[k]) / (p4m - 1);
					}
					rmax--;
				}
				res[u] = Daprox[0];
			}
			u++;
		}
    }

	return res;
}


int hessian(TFunc func,
			const ValArr& x,
			ValMat &H,
			const std::string& method,
			const TArgs& method_args)
{

    // if(method=="complex"){ # Complex step hessian
    //   args <- list(eps=1e-4, d=0.1,
    //      zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
    //   args[names(method.args)] <- method.args
    //   # the CSD part of this uses eps=.Machine$double.eps
    //   # but the jacobian is Richardson and uses method.args
    //   fn <- function(x, ...){
    //           grad(func=func, x=x, method="complex", side=NULL,
    //	   method.args=list(eps=.Machine$double.eps), ...)
    //	   }
    //
    //   return(jacobian(func=fn, x=x, method="Richardson", side=NULL,
    //                    method.args=args, ...))
    //   }
    // else
    int nx = func.DimX();
	int ny = func.DimY();
	
    if(method != "Richardson")
		stop("method not implemented.");
    if(1 != ny)
		stop("Richardson method for hessian assumes a scalar valued function.");

    // H
    H.resize(nx, nx);
	H.setZero();

    ValArr D = genD(func, x, method, method_args);

    // H <- diag(NA,length(x))
    int u = nx;
    for(int i = 0; i < nx; i++) {
		for(int j = 0; j <= i; j++) {
			H(i, j) = D[u++];
		}
    }
    H = H + H.transpose().eval();

    // diag(H) <- diag(H)/2
    for(int i = 0; i < nx; i++) {
		H(i, i) = H(i, i) / 2.0;
    }

    return 1;
}

ValMat hessian(TFunc func,
    const ValArr& x,
    const std::string& method,
    const TArgs& method_args)
{
	int nx = func.DimX();
	ValMat H(nx, nx);
	hessian(func, x, H, method, method_args);
	return H;
}


template<class T>
bool doTest(const T &res, const T &ref, double tol, const std::string &tname, bool &testSuite){
	double err = abs(res - ref).maxCoeff();
	bool testOK = err < tol;
	if(!testOK){
		std::cout << tname << " FAILED: err=" << err << " > " << tol << std::endl;
		testSuite = false;
	} else {
		std::cout << tname << " PASSED: err=" << err << " < " << tol << std::endl;
	}
	return testOK;
}
	
bool testDeriv() {
	bool testSuite(true);
	//bool testOK(true);
	double My_PI(3.14159265358979323846);
	//double err;
	int Ny(1);
	int Nx(11);
    //
	ValArr side_par(0);
	ValArr x_par(Nx);
    TArgs args = std::make_tuple(0.0001, 0.0001, 1.781e-5, 4, 2, false);

    ValArr res(Nx);
	ValArr ref(Nx);
	
//	[1]  1.000000  0.809017  0.309017 -0.309017 -0.809017 -1.000000 -0.809017 -0.309017  0.309017  0.809017  1.000000
    x_par << 0.0, 0.6283185, 1.2566371, 1.8849556, 2.5132741, 3.1415927, 3.7699112, 4.3982297, 5.0265482, 5.6548668, 6.2831853;
    TFunc ff(Nx, Ny, [](const ValArr &x_dom){ return sin(x_dom).sum(); });
	ref = cos(x_par);

    std::cout << std::endl << "Gradient: Simple, f=sum(sin(x))";
	res = gradient(ff, x_par, "simple", side_par, args);
    std::cout << std::endl << res << std::endl;
	doTest<ValArr>(res, ref, 1e-4, "grad01 test 4", testSuite);
	
    std::cout << std::endl << "Gradient: Richardson, f=sum(sin(x))";
    res = gradient(ff, x_par, "Richardson", side_par, args);
    std::cout << std::endl << res << std::endl;
	doTest<ValArr>(res, ref, 1e-10, "grad01 test 3", testSuite);
	
	Nx = 3;
	x_par.resize(Nx);
	ValMat resH(Nx,Nx);
	ValMat refH(Nx,Nx);
//	     [,1] [,2] [,3]
//	[1,]    0    0    0
//	[2,]    0    0    0
//	[3,]    0    0    0
	std::cout << std::endl << "Hessian: f=sum(sin(x))";
	x_par << 0.0, 3.141593, 6.283185;
	ff.SetDimX(Nx);
	refH.setZero();
	resH = hessian(ff, x_par, "Richardson", args);
	std::cout << std::endl << resH << std::endl;
	doTest<ValMat>(resH, refH, 1e-10, "hessian test 2", testSuite);
	
//			[,1]     [,2]     [,3]
//	[1,] 29.55622    0.000     0.00
//	[2,]  0.00000 1613.715     0.00
//	[3,]  0.00000    0.000 88105.86
	std::cout << std::endl << "Hessian: f=sum(exp(2x))";
	x_par << 1.0, 3.0, 5.0;
	TFunc ff4(Nx, Ny, [](const ValArr &x_dom){ return exp(2*x_dom).sum(); });
	std::get<_D>(args) = 0.01;
	ValArr v1 = 4*exp(2*x_par);
	refH(0,0) = v1(0);
	refH(1,1) = v1(1);
	refH(2,2) = v1(2);
	resH = hessian(ff4, x_par, "Richardson", args);
	std::cout << std::endl << resH << std::endl;
	doTest<ValMat>(resH, refH, 1e-5, "hessian test 7", testSuite);
	
	Nx = 1;
	x_par.resize(Nx);
	resH.resize(Nx, Nx);
	refH.resize(Nx, Nx);
	std::get<_D>(args) = 0.0001;
	
//			[,1]
//	[1,] -0.7071068
	std::cout << std::endl << "Hessian: f=sin(x)";
	x_par << 0.25*My_PI;
	TFunc ff2(Nx, Ny, [](const ValArr &x_dom){ return sin(x_dom[0]); });
	std::get<_D>(args) = 0.0001;
	refH(0,0) = sin(x_par(0)+My_PI);
	resH = hessian(ff2, x_par, "Richardson", args);
	std::cout << std::endl << resH << std::endl;
	doTest<ValMat>(resH, refH, 1e-4, "hessian test 3", testSuite);
	
//			[,1]
//	[1,] 29.55622
	std::cout << std::endl << "Hessian: f=exp(2x)";
	x_par << 1.0;
	TFunc ff3(Nx, Ny, [](const ValArr &x_dom){ return exp(2*(x_dom[0])); });
	refH(0,0) = 4*exp(2*x_par(0));
	resH = hessian(ff3, x_par, "Richardson", args);
	std::cout << std::endl << resH << std::endl;
	doTest<ValMat>(resH, refH, 1e-3, "hessian test 5", testSuite);
	
	return testSuite;
}

#endif