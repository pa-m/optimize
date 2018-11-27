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
