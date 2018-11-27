package optimize

import (
	"log"
	"math"
)

// PowellMinimizer minimizes a scalar function of multidimensionnal x using modified Powell algorithm
// (see fmin_powell in scipy.optimize)
type PowellMinimizer struct {
	Callback        func([]float64)
	Xtol, Ftol      float64
	MaxIter, MaxFev int
	Logger          *log.Logger
}

// NewPowellMinimizer return a PowellMinimizer with default tolerances
func NewPowellMinimizer() (pm *PowellMinimizer) {
	pm = &PowellMinimizer{Xtol: 1e-4, Ftol: 1e-4}
	return
}

// Minimize minimizes f starting at x0
func (pm *PowellMinimizer) Minimize(f func([]float64) float64, x0 []float64) {
	minimizePowell(f, x0, pm.Callback, pm.Xtol, pm.Ftol, pm.MaxIter, pm.MaxFev, pm.Logger)
}

// Minimization of scalar function of one or more variables using the
// modified Powell algorithm.
// Options
// -------
// disp : bool
//     Set to True to print convergence messages.
// xtol : float
//     Relative error in solution `xopt` acceptable for convergence.
// ftol : float
//     Relative error in ``fun(xopt)`` acceptable for convergence.
// maxiter, maxfev : int
//     Maximum allowed number of iterations and function evaluations.
//     Will default to ``N*1000``, where ``N`` is the number of
//     variables, if neither `maxiter` or `maxfev` is set. If both
//     `maxiter` and `maxfev` are set, minimization will stop at the
//     first reached.
// direc : ndarray
//     Initial set of direction vectors for the Powell method.
func minimizePowell(
	f func([]float64) float64, x0 []float64, callback func([]float64),
	xtol, ftol float64,
	maxiter, maxfev int,
	disp *log.Logger) []float64 {
	type float = float64
	var (
		fval, fx, delta, fx2, bnd, t, temp float
		x1, x2, direc, direc1              []float
		bigind, warnflag                   int
	)
	abs := func(x float) float {
		if x < 0 {
			return -x
		}
		return x
	}
	// # we need to use a mutable object here that we can update in the
	// # wrapper function
	fcalls := 0
	fun := func(x []float) float {
		y := f(x)
		fcalls++
		return y
	}
	x := make([]float64, len(x0))
	N := len(x)
	//# If neither are set, then set both to default
	if maxiter <= 0 && maxfev <= 0 {
		maxiter = N * 1000
		maxfev = N * 1000
	} else if maxiter <= 0 {
		// # Convert remaining Nones, to np.inf, unless the other is np.inf, in
		// # which case use the default to avoid unbounded iteration
		if maxfev == math.MaxInt64 {
			maxiter = N * 1000
		} else {
			maxiter = math.MaxInt64
		}
	} else if maxfev <= 0 {
		if maxiter == math.MaxInt64 {
			maxfev = N * 1000
		} else {
			maxfev = math.MaxInt64
		}
	}
	// direc is used as a matrix direc[i,j]:=direc[i*N+j]
	direc = make([]float, N*N)
	direc1 = make([]float, N)
	for i := 0; i < N; i++ {
		direc[i*N+i] = 1
	}

	fval = fun(x)
	x1, x2 = make([]float64, N), make([]float64, N)
	copy(x1, x)
	iter := 0
	ilist := make([]int, N)
	for i := range ilist {
		ilist[i] = i
	}
	for {
		fx = fval
		bigind = 0
		delta = 0.0
		for _, i := range ilist {
			direc1 = direc[i*N : i*N+N]
			fx2 = fval
			fval, x, direc1 = linesearchPowell(fun, x, direc1, xtol*100)
			if (fx2 - fval) > delta {
				delta = fx2 - fval
				bigind = i
			}
		}
		iter++
		if callback != nil {
			callback(x)
		}
		bnd = ftol*(abs(fx)+abs(fval)) + 1e-20
		if 2.0*(fx-fval) <= bnd {
			break
		}
		if fcalls >= maxfev {
			break
		}
		if iter >= maxiter {
			break
		}
		//# Construct the extrapolated point
		// direc1 = x - x1
		// x2 = 2*x - x1
		// x1 = x.copy()
		for i, xi := range x {
			direc1[i] = xi - x1[i]
			x2[i] = 2*xi - x1[i]
			x1[i] = xi
		}
		fx2 = fun(x2)

		if fx > fx2 {
			t = 2.0 * (fx + fx2 - 2.0*fval)
			temp = (fx - fval - delta)
			t *= temp * temp
			temp = fx - fx2
			t -= delta * temp * temp
			if t < 0.0 {
				fval, x, direc1 = linesearchPowell(fun, x, direc1, xtol*100)
				//direc[bigind] = direc[-1]
				copy(direc[bigind*N:bigind*N+N], direc[(N-1)*N:N*N])
				//direc[-1] = direc1
				copy(direc[(N-1)*N:N*N], direc1)
			}
		}

	}
	warnflag = 0
	if fcalls >= maxfev {
		warnflag = 1
		//msg = _status_message['maxfev']
		msg := "maxfev"
		if disp != nil {
			disp.Println("Warning: " + msg)
		}
	} else if iter >= maxiter {
		warnflag = 2
		//msg = _status_message['maxiter']
		msg := "maxiter"
		if disp != nil {
			disp.Println("Warning: " + msg)
		}
	} else {
		//msg = _status_message['success']
		if disp != nil {
			disp.Printf("Success. Current function value: %.7g Iterations: %d Function evaluations: %d", fval, iter, fcalls)
		}
	}
	//x = squeeze(x)
	_ = warnflag
	return x
}

// Line-search algorithm using fminbound. Find the minimum of the function ``func(x0+ alpha*direc)``.
func linesearchPowell(
	fun func([]float64) float64,
	p, xi []float64,
	tol float64) (float64, []float64, []float64) {
	type float = float64
	myfunc := func(alpha float) float {

		//return fun(p + alpha*xi)
		xtmp := make([]float, len(p))
		for i, p1 := range p {
			xtmp[i] = p1 + alpha*xi[i]
		}
		return fun(xtmp)
	}
	b := newBrentMiner(myfunc, tol, 500)
	b.optimize()
	alphaMin, fret := b.Xmin, b.Fval
	//xi = alpha_min*xi
	//return squeeze(fret), p + xi, xi
	pPlusXi := make([]float, len(p))
	for i := range p {
		xi[i] *= alphaMin
		pPlusXi[i] = p[i] + xi[i]
	}

	return fret, pPlusXi, xi
}

type bracketer struct {
	growLimit float64
	maxIter   int
}

// Bracket the minimum of the function.
// Given a function and distinct initial points, search in the
// downhill direction (as defined by the initital points) and return
// new points xa, xb, xc that bracket the minimum of the function
// f(xa) > f(xb) < f(xc). It doesn't always mean that obtained
// solution will satisfy xa<=x<=xb
func (b bracketer) bracket(f func(float64) float64, xa0, xb0 float64) (xa, xb, xc, fa, fb, fc float64, funcalls int) {
	var (
		tmp1, tmp2, val, denom, w, wlim, fw float64
		iter                                int
	)
	_gold := 1.618034 //# golden ratio: (1.0+sqrt(5.0))/2.0
	_verysmallNum := 1e-21
	xa, xb = xa0, xb0
	fa, fb = f(xa), f(xb)
	if fa < fb {
		xa, xb = xb, xa
		fa, fb = fb, fa
	}
	xc = xb + _gold*(xb-xa)
	fc = f(xc)
	funcalls = 3
	iter = 0
	for fc < fb {
		tmp1 = (xb - xa) * (fb - fc)
		tmp2 = (xb - xc) * (fb - fa)
		val = tmp2 - tmp1
		if math.Abs(val) < _verysmallNum {
			denom = 2.0 * _verysmallNum
		} else {
			denom = 2.0 * val
		}
		w = xb - ((xb-xc)*tmp2-(xb-xa)*tmp1)/denom
		wlim = xb + b.growLimit*(xc-xb)
		if iter > b.maxIter {
			log.Fatal("Too many iterations.")
		}
		iter++
		if (w-xc)*(xb-w) > 0.0 {
			fw = f(w)
			funcalls++
			if fw < fc {
				xa = xb
				xb = w
				fa = fb
				fb = fw
				return xa, xb, xc, fa, fb, fc, funcalls
			} else if fw > fb {
				xc = w
				fc = fw
				return xa, xb, xc, fa, fb, fc, funcalls
			}
			w = xc + _gold*(xc-xb)
			fw = f(w)
			funcalls++
		} else if (w-wlim)*(wlim-xc) >= 0.0 {
			w = wlim
			fw = f(w)
			funcalls++
		} else if (w-wlim)*(xc-w) > 0.0 {
			fw = f(w)
			funcalls++
			if fw < fc {
				xb = xc
				xc = w
				w = xc + _gold*(xc-xb)
				fb = fc
				fc = fw
				fw = f(w)
				funcalls++
			}
		} else {
			w = xc + _gold*(xc-xb)
			fw = f(w)
			funcalls++
		}
		xa = xb
		xb = xc
		xc = w
		fa = fb
		fb = fc
		fc = fw
	}
	return xa, xb, xc, fa, fb, fc, funcalls
}

// brentMinimizer is the translation of class Brent in scipy/optimize/optimize.py
type brentMinimizer struct {
	Func           func(float64) float64
	Tol            float64
	Maxiter        int
	mintol         float64
	cg             float64
	Xmin           float64
	Fval           float64
	Iter, Funcalls int
	brack          []float64
	bracketer
}

func newBrentMiner(fun func(float64) float64, tol float64, maxiter int) *brentMinimizer {
	return &brentMinimizer{
		Func:      fun,
		Tol:       tol,
		Maxiter:   maxiter,
		mintol:    1.0e-11,
		cg:        0.3819660,
		bracketer: bracketer{growLimit: 110, maxIter: 1000},
	}
}
func (bm *brentMinimizer) setBracket(brack []float64) {
	bm.brack = make([]float64, len(brack))
	copy(bm.brack, brack)
}
func (bm *brentMinimizer) getBracketInfo() (float64, float64, float64, float64, float64, float64, int) {
	fun := bm.Func
	brack := bm.brack
	var xa, xb, xc float64
	var fa, fb, fc float64
	var funcalls int
	switch len(brack) {
	case 0:
		xa, xb, xc, fa, fb, fc, funcalls = bm.bracketer.bracket(fun, 0, 1)
	case 2:
		xa, xb, xc, fa, fb, fc, funcalls = bm.bracketer.bracket(fun, brack[0], brack[1])
	case 3:
		xa, xb, xc = brack[0], brack[1], brack[2]
		if xa > xc {
			xa, xc = xc, xa
		}
		fa, fb, fc = fun(xa), fun(xb), fun(xc)
		if !((fb < fa) && (fb < fc)) {
			log.Fatal("not a brackeding interval")
		}
		funcalls = 3
	}
	return xa, xb, xc, fa, fb, fc, funcalls
}

func (bm *brentMinimizer) optimize() {
	var (
		xa, xb, xc, fb, _mintol, _cg, x, fx, v, fv, w, fw, a, b, deltax, tol1, tol2, xmid, rat, tmp1, tmp2, p, dxTemp, u, fu float64
		funcalls, iter                                                                                                       int
	)
	//# set up for optimization
	f := bm.Func
	xa, xb, xc, _, fb, _, funcalls = bm.getBracketInfo()
	_mintol = bm.mintol
	_cg = bm.cg
	// #################################
	// #BEGIN CORE ALGORITHM
	//#################################
	//x = w = v = xb
	v, w, x = xb, xb, xb
	//fw = fv = fx = func(*((x,) + self.args))
	fx = fb
	fv, fw = fx, fx
	if xa < xc {
		a = xa
		b = xc
	} else {
		a = xc
		b = xa
	}
	deltax = 0.0
	funcalls++
	iter = 0
	for iter < bm.Maxiter {
		tol1 = bm.Tol*math.Abs(x) + _mintol
		tol2 = 2.0 * tol1
		xmid = 0.5 * (a + b)
		//# check for convergence
		if math.Abs(x-xmid) < (tol2 - 0.5*(b-a)) {
			break
		}
		// # XXX In the first iteration, rat is only bound in the true case
		// # of this conditional. This used to cause an UnboundLocalError
		// # (gh-4140). It should be set before the if (but to what?).
		if math.Abs(deltax) <= tol1 {
			if x >= xmid {
				deltax = a - x //# do a golden section step
			} else {
				deltax = b - x
			}
			rat = _cg * deltax
		} else { //# do a parabolic step
			tmp1 = (x - w) * (fx - fv)
			tmp2 = (x - v) * (fx - fw)
			p = (x-v)*tmp2 - (x-w)*tmp1
			tmp2 = 2.0 * (tmp2 - tmp1)
			if tmp2 > 0.0 {
				p = -p
			}
			tmp2 = math.Abs(tmp2)
			dxTemp = deltax
			deltax = rat
			//# check parabolic fit
			if (p > tmp2*(a-x)) && (p < tmp2*(b-x)) &&
				(math.Abs(p) < math.Abs(0.5*tmp2*dxTemp)) {
				rat = p * 1.0 / tmp2 //# if parabolic step is useful.
				u = x + rat
				if (u-a) < tol2 || (b-u) < tol2 {
					if xmid-x >= 0 {
						rat = tol1
					} else {
						rat = -tol1
					}
				}
			} else {
				if x >= xmid {
					deltax = a - x //# if it's not do a golden section step
				} else {
					deltax = b - x
				}
				rat = _cg * deltax
			}
		}
		if math.Abs(rat) < tol1 { //# update by at least tol1
			if rat >= 0 {
				u = x + tol1
			} else {
				u = x - tol1
			}
		} else {
			u = x + rat
		}
		fu = f(u) //# calculate new output value
		funcalls++

		if fu > fx { //# if it's bigger than current
			if u < x {
				a = u
			} else {
				b = u
			}
			if (fu <= fw) || (w == x) {
				v = w
				w = u
				fv = fw
				fw = fu
			} else if (fu <= fv) || (v == x) || (v == w) {
				v = u
				fv = fu
			}
		} else {
			if u >= x {
				a = x
			} else {
				b = x
			}
			v = w
			w = x
			x = u
			fv = fw
			fw = fx
			fx = fu
		}
		iter++
	}
	// #################################
	// #END CORE ALGORITHM
	// #################################
	bm.Xmin, bm.Fval, bm.Iter, bm.Funcalls = x, fx, iter, funcalls
}
