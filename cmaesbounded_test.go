package optimize

// func ExampleCmaEsCholB() {
// 	problem := optimize.Problem{
// 		Func: func(x []float64) float64 {
// 			return x[0]*x[0] + x[1]*x[1]
// 		},
// 	}
// 	initX := []float64{1, 1}
// 	method := &CmaEsCholB{Xmin: []float64{.1, math.Inf(-1)}}
// 	method.Src = rand.NewSource(uint64(1))
// 	settings := &optimize.Settings{FuncEvaluations: 500}

// 	res, err := optimize.Minimize(problem, initX, settings, method)
// 	if err != nil {
// 		panic(err)
// 	}
// 	//fmt.Printf("%#v\n", res)
// 	if math.Abs(res.Location.X[0]-.1) > 1e-2 || math.Abs(res.Location.X[1]-.0) > 1e-2 {
// 		fmt.Printf("%.5f", res.Location.X)

// 	}
// 	// Output:
// }
