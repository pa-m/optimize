package optimize

import "gonum.org/v1/gonum/mat"

const (
	nonpositiveDimension string = "optimize: non-positive input dimension"
	negativeTasks        string = "optimize: negative input number of tasks"
)

// resize takes x and returns a slice of length dim. It returns a resliced x
// if cap(x) >= dim, and a new slice otherwise.
func resize(x []float64, dim int) []float64 {
	if dim > cap(x) {
		return make([]float64, dim)
	}
	return x[:dim]
}

func resizeSymDense(m *mat.SymDense, dim int) *mat.SymDense {
	if m == nil || cap(m.RawSymmetric().Data) < dim*dim {
		return mat.NewSymDense(dim, nil)
	}
	return mat.NewSymDense(dim, m.RawSymmetric().Data[:dim*dim])
}
