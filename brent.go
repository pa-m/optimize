package main

import (
	"errors"
	"fmt"
	"log"
	"math"
)

// Brent find zero of f using Brent's method
// see https://en.wikipedia.org/wiki/Brent%27s_method
// logger may be nil
func Brent(a, b, tol float64, f func(float64) float64, logger *log.Logger) (float64, error) {
	type float = float64
	abs := func(x float) float {
		if x < 0 {
			return -x
		}
		return x
	}
	it := 0
	// Entrée de a, b, et un pointeur vers une sous-routine pour f
	// calculer f(a)
	// calculer f(b)
	fa, fb := f(a), f(b)
	// si f(a) f(b) >= 0 alors sortie (erreur) fin si
	if fa*fb >= 0 {
		return math.NaN(), errors.New("brent: f(a) f(b) >= 0")
	}
	// si |f(a)| < |f(b)| alors échanger (a,b) fin si
	if abs(fa) < abs(fb) {
		a, fa, b, fb = b, fb, a, fa
	}
	// c := a
	c, fc := a, fa
	var d, s, fs float
	// mflag := vrai
	mflag := true
	// répéter jusqu'à ce que f(b) = 0 ou |b − a| soit suffisamment petit (convergence)
	for fb != 0 && abs(b-a) > tol {
		if logger != nil {
			logger.Printf("%d (a%d,f(a%d))=(%.5g, %.5g) and  (b%d,f(b%d))=%.5g,%.5g ", it+1, it, it, a, fa, it, it, b, fb)
		}
		it++
		if it == 1000 {
			return math.NaN(), fmt.Errorf("brent: it=%d", it)
		}
		//     si f(a) ≠ f(c) et f(b) ≠ f(c) alors
		//         s := a f ( b ) f ( c ) ( f ( a ) − f ( b ) ) ( f ( a ) − f ( c ) ) + b f ( a ) f ( c ) ( f ( b ) − f ( a ) ) ( f ( b ) − f ( c ) ) + c f ( a ) f ( b ) ( f ( c ) − f ( a ) ) ( f ( c ) − f ( b ) ) {\displaystyle s:={\frac {af(b)f(c)}{(f(a)-f(b))(f(a)-f(c))}}+{\frac {bf(a)f(c)}{(f(b)-f(a))(f(b)-f(c))}}+{\frac {cf(a)f(b)}{(f(c)-f(a))(f(c)-f(b))}}} s:={\frac {af(b)f(c)}{(f(a)-f(b))(f(a)-f(c))}}+{\frac {bf(a)f(c)}{(f(b)-f(a))(f(b)-f(c))}}+{\frac {cf(a)f(b)}{(f(c)-f(a))(f(c)-f(b))}} (interpolation quadratique inverse)
		//     sinon
		//         s := b − f ( b ) b − a f ( b ) − f ( a ) {\displaystyle s:=b-f(b){\frac {b-a}{f(b)-f(a)}}} s:=b-f(b){\frac {b-a}{f(b)-f(a)}} (règle de la sécante)
		//     fin si
		if fa != fc && fb != fc {
			s = a*fb*fc/(fa-fb)/(fa-fc) +
				b*fa*fc/(fb-fa)/(fb-fc) +
				c*fa*fb/(fc-fa)/(fc-fb)
		} else {
			s = b - fb*(b-a)/(fb-fa)
		}

		//     si s n'est pas entre (3a + b)/4 et b ou (mflag est vrai et |s−b| ≥ |b−c| / 2) ou (mflag est faux et |s−b| ≥ |c−d| / 2) alors
		//         s := a + b 2 {\displaystyle s:={\frac {a+b}{2}}} s:={\frac {a+b}{2}}
		//         mflag := vrai
		//     sinon
		//         mflag := faux
		//     fin si
		between := ((3*a+b)/4 <= s && s <= b) || ((3*a+b)/4 >= s && s >= b)
		var ineq bool
		if between {
			if mflag {
				ineq = abs(s-b) < abs(b-c)/2
			} else {
				ineq = abs(s-b) < abs(c-d)/2
			}
		}

		if (!between) || !ineq {
			s = (a + b) / 2
			mflag = true
		} else {
			mflag = false
		}

		//     calculer f(s)
		fs = f(s)
		//     d := c
		//     c := b
		d = c
		c, fc = b, fb
		//     si f(a) f(s) < 0 alors b := s sinon a := s fin si
		if fa*fs < 0 {
			b, fb = s, fs
		} else {
			a, fa = s, fs
		}
		//     si |f(a)| < |f(b)| alors échange (a,b) fin si
		if abs(fa) < abs(fb) {
			a, fa, b, fb = b, fb, a, fa
		}
		// fin répéte
	}
	if logger != nil {
		logger.Printf("%d (a%d,f(a%d))=(%.5g, %.5g) and  (b%d,f(b%d))=%.5g,%.5g ", it+1, it, it, a, fa, it, it, b, fb)
	}
	// sortir b (renvoie de la racine)
	return b, nil
}

// Bissection find zero of f using Bissection's method
// logger may be nil
func Bissection(a, b, tol float64, f func(float64) float64, logger *log.Logger) (float64, error) {
	type float = float64
	abs, NaN := math.Abs, math.NaN()
	it := 0
	// Entrée de a, b, et un pointeur vers une sous-routine pour f
	// calculer f(a)
	// calculer f(b)
	fa, fb := f(a), f(b)
	// si f(a) f(b) >= 0 alors sortie (erreur) fin si
	if fa*fb >= 0 {
		return NaN, errors.New("brent: f(a) f(b) >= 0")
	}
	// si |f(a)| < |f(b)| alors échanger (a,b) fin si
	if abs(fa) < abs(fb) {
		a, fa, b, fb = b, fb, a, fa
	}
	var s, fs float
	// répéter jusqu'à ce que f(b) = 0 ou |b − a| soit suffisamment petit (convergence)
	for fb != 0 && abs(b-a) > tol {
		if logger != nil {
			logger.Printf("%d a,fa=%.5g, %.5g b,fb=%.5g,%.5g\n", it, a, fa, b, fb)
		}
		it++
		s = (a + b) / 2
		fs = f(s)
		//     si f(a) f(s) < 0 alors b := s sinon a := s fin si
		if fa*fs < 0 {
			b, fb = s, fs
		} else {
			a, fa = s, fs
		}
		//     si |f(a)| < |f(b)| alors échange (a,b) fin si
		if abs(fa) < abs(fb) {
			a, fa, b, fb = b, fb, a, fa
		}
		// fin répéte
	}
	if logger != nil {
		logger.Printf("%d a,fa=%.5g, %.5g b,fb=%.5g,%.5g\n", it, a, fa, b, fb)
	}
	// sortir b (renvoie de la racine)
	return b, nil
}
