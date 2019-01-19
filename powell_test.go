package optimize

import (
	"fmt"
	"log"
	"math"
	"os"
)

func ExamplePowellMinimizer() {
	pm := NewPowellMinimizer()
	pm.Callback = func(x []float64) {
		fmt.Printf("%.5f\n", x)
	}
	pm.Logger = log.New(os.Stdout, "", 0)

	pm.Minimize(
		func(x []float64) float64 { return -math.Exp(1 / (1 + x[0]*x[0] + x[1]*x[1])) },
		[]float64{10, 20},
	)
	// Output:
	// [-0.02748 -0.02037]
	// [0.00818 -0.00407]
	// [0.00154 -0.00337]
	// Success. Current function value: -2.718245 Iterations: 3 Function evaluations: 69
}
