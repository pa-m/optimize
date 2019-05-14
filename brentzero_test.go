package optimize

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
	_, err := Brent(a, b, 1e-9, f, log.New(os.Stdout, "", 0))
	if err != nil {
		fmt.Println(err.Error())
	}

	// // Output:
	// 1 (a0,f(a0))=(-4, -25) and  (b0,f(b0))=1.3333,0.48148
	// 2 (a1,f(a1))=(-4, -25) and  (b1,f(b1))=1.2326,0.22891
	// 3 (a2,f(a2))=(-4, -25) and  (b2,f(b2))=1.1421,0.083582
	// 4 (a3,f(a3))=(-4, -25) and  (b3,f(b3))=-1.429,9.2689
	// 5 (a4,f(a4))=(-4, -25) and  (b4,f(b4))=-2.7145,3.9393
	// 6 (a5,f(a5))=(-3.3572, -6.7825) and  (b5,f(b5))=-2.7145,3.9393
	// 7 (a6,f(a6))=(-2.7145, 3.9393) and  (b6,f(b6))=-3.0359,-0.58418
	// 8 (a7,f(a7))=(-3.0359, -0.58418) and  (b7,f(b7))=-2.9944,0.089961
	// 9 (a8,f(a8))=(-3.0359, -0.58418) and  (b8,f(b8))=-2.9999,0.0015995
	// 10 (a9,f(a9))=(-2.9999, 0.0015995) and  (b9,f(b9))=-3,-1.3912e-07
	// 11 (a10,f(a10))=(-3, -1.3912e-07) and  (b10,f(b10))=-3,6.9562e-12
	// 12 (a11,f(a11))=(-3, 6.9562e-12) and  (b11,f(b11))=-3,0
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
	_, err := Bissection(a, b, 1e-9, f, log.New(os.Stdout, "", 0))
	if err != nil {
		panic(err)
	}
	// Output:
	// 0 a,fa=-4, -25 b,fb=1.3333,0.48148
	// 1 a,fa=-4, -25 b,fb=-1.3333,9.0741
	// 2 a,fa=-4, -25 b,fb=-2.6667,4.4815
	// 3 a,fa=-3.3333, -6.2593 b,fb=-2.6667,4.4815
	// 4 a,fa=-2.6667, 4.4815 b,fb=-3,0
}
