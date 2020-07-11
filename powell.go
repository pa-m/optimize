package optimize

import (
	"log"
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
	const MaxInt = (int)(^uint(0) >> 1)
	//# If neither are set, then set both to default
	N := len(x0)
	if pm.MaxIter <= 0 && pm.MaxFev <= 0 {
		pm.MaxIter = N * 1000
		pm.MaxFev = N * 1000
	} else if pm.MaxIter <= 0 {
		// # Convert remaining Nones, to np.inf, unless the other is np.inf, in
		// # which case use the default to avoid unbounded iteration
		if pm.MaxFev == MaxInt {
			pm.MaxIter = N * 1000
		} else {
			pm.MaxIter = MaxInt
		}
	} else if pm.MaxFev <= 0 {
		if pm.MaxIter == MaxInt {
			pm.MaxFev = N * 1000
		} else {
			pm.MaxFev = MaxInt
		}
	}
	fnMaxIter := func(iter int) bool { return iter >= pm.MaxIter }
	fnMaxFev := func(fcalls int) bool { return fcalls >= pm.MaxFev }
	minimizePowell(f, x0, pm.Callback, pm.Xtol, pm.Ftol, fnMaxIter, fnMaxFev, pm.Logger)
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
	f func([]float64) float64,
	x0 []float64,
	callback func([]float64),
	xtol, ftol float64,
	fnMaxIter func(int) bool, fnMaxFev func(int) bool,
	disp *log.Logger) ([]float64, int) {
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
	if fnMaxIter == nil {
		fnMaxIter = func(int) bool { return false }
	}
	if fnMaxFev == nil {
		fnMaxIter = func(int) bool { return false }
	}
	// # we need to use a mutable object here that we can update in the
	// # wrapper function
	fcalls := 0
	fun := func(x []float) float {
		y := f(x)
		fcalls++
		return y
	}
	fnMaxFevSub := func(funcalls int) bool { return fnMaxFev(fcalls + funcalls) }
	if callback == nil {
		callback = func(x []float64) {}
	}
	N := len(x0)
	x := make([]float64, N)
	copy(x, x0)

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
			fval, x, direc1 = linesearchPowell(fun, x, direc1, xtol*100, fnMaxFevSub)
			if (fx2 - fval) > delta {
				delta = fx2 - fval
				bigind = i
			}
		}
		iter++
		callback(x)
		bnd = ftol*(abs(fx)+abs(fval)) + 1e-20
		if 2.0*(fx-fval) <= bnd {
			break
		}
		if fnMaxFev(fcalls) {
			break
		}
		if fnMaxIter(iter) {
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
				fval, x, direc1 = linesearchPowell(fun, x, direc1, xtol*100, fnMaxFevSub)
				//direc[bigind] = direc[-1]
				copy(direc[bigind*N:bigind*N+N], direc[(N-1)*N:N*N])
				//direc[-1] = direc1
				copy(direc[(N-1)*N:N*N], direc1)
			}
		}

	}
	warnflag = 0
	if fnMaxFev(fcalls) {
		// FunctionEvaluationLimit
		warnflag = 1
		//msg = _status_message['maxfev']
		msg := "maxfev"
		if disp != nil {
			disp.Println("Warning: " + msg)
		}
	} else if fnMaxIter(iter) {
		// IterationLimit
		warnflag = 2
		//msg = _status_message['maxiter']
		msg := "maxiter"
		if disp != nil {
			disp.Println("Warning: " + msg)
		}
	} else {
		// Success,MethodConverge ?
		//msg = _status_message['success']
		if disp != nil {
			disp.Printf("Success. Current function value: %.7g Iterations: %d Function evaluations: %d", fval, iter, fcalls)
		}
	}
	return x, warnflag
}

// Line-search algorithm using fminbound. Find the minimum of the function ``func(x0+ alpha*direc)``.
func linesearchPowell(
	fun func([]float64) float64,
	p, xi []float64,
	tol float64,
	fnMaxFev func(int) bool,
) (float64, []float64, []float64) {
	type float = float64
	myfunc := func(alpha float) float {

		//return fun(p + alpha*xi)
		xtmp := make([]float, len(p))
		for i, p1 := range p {
			xtmp[i] = p1 + alpha*xi[i]
		}
		return fun(xtmp)
	}

	alphaMin, fret, _, _ := NewBrentMinimizer(myfunc, tol, 500, fnMaxFev).Optimize()
	//xi = alpha_min*xi
	//return squeeze(fret), p + xi, xi
	pPlusXi := make([]float, len(p))
	for i := range p {
		xi[i] *= alphaMin
		pPlusXi[i] = p[i] + xi[i]
	}

	return fret, pPlusXi, xi
}
