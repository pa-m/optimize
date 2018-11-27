package optimize

import (
	"fmt"
	"log"
	"math"
	"os"
)

func ExamplePowellMinimizer() {

	square := func(x float64) float64 {
		return x * x
	}
	pm := NewPowellMinimizer()
	pm.Callback = func(x []float64) { fmt.Printf("%.5g\n", x) }
	pm.Logger = log.New(os.Stdout, "", 0)

	pm.Minimize(
		func(x []float64) float64 { return math.Log(square(x[0]-2) + square(x[1]-3) + 4) },
		[]float64{10, 20},
	)

	// Output:
	// [1.9998 3.0111]
	// [1.9949 3.0034]
	// Success. Current function value: 1.386304 Iterations: 2 Function evaluations: 40

}
