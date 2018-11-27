package main

import (
	"fmt"
	"log"
	"os"
)

func ExampleBrent() {
	//On cherche à identifier une racine de f(x) = (x + 3)(x − 1)2
	f := func(x float64) float64 {
		xless1 := x - 1
		y := (x + 3) * xless1 * xless1
		return y
	}
	//On prend [a0; b0] = [−4; 4/3]
	a, b := -4.0, 4./3.
	xstar, err := Brent(a, b, 1e-9, f, log.New(os.Stdout, "", 0))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("%.5g", xstar)
	// Output:
	// 0 a,fa=       -4,       -25 b,fb=  1.33333, 0.481481
	// 1 a,fa=       -4,       -25 b,fb=  1.23256, 0.228911
	// 2 a,fa=       -4,       -25 b,fb= -1.38372,   9.1839
	// 3 a,fa=       -4,       -25 b,fb= -2.69186,  4.19989
	// 4 a,fa= -3.34593,  -6.53362 b,fb= -2.69186,  4.19989
	// 5 a,fa= -3.34593,  -6.53362 b,fb= -2.94779, 0.813697
	// 6 a,fa= -2.94779,  0.813697 b,fb= -3.00248,-0.0397038
	// 7 a,fa= -3.00248, -0.0397038 b,fb= -2.99993,0.00105477
	// 8 a,fa= -3.00248, -0.0397038 b,fb=       -3,1.30594e-06
	// 9 a,fa=       -3, 1.30594e-06 b,fb=       -3,-9.9476e-14
	// 10 a,fa=-3, -9.9476e-14 b,fb=-3,0
	// -3

}

func ExampleBissection() {
	//On cherche à identifier une racine de f(x) = (x + 3)(x − 1)2
	f := func(x float64) float64 {
		xless1 := x - 1
		y := (x + 3) * xless1 * xless1
		return y
	}
	//On prend [a0; b0] = [−4; 4/3]
	a, b := -4.0, 4./3.
	xstar, err := Bissection(a, b, 1e-9, f, log.New(os.Stdout, "", 0))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("%.5g", xstar)
	// Output:
	// -3
}
