// adapted from gonum.org/v1/gonum/optimize/cmaes.go. added ensureBounds

// Copyright ©2017 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
/* BSD license for code copied from gonum/optimize/cmaes.go (all except sendTask,ensureBounds)
Copyright ©2013 The Gonum Authors. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the gonum project nor the names of its authors and
      contributors may be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

package optimize

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/optimize"

	"golang.org/x/exp/rand"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat/distmv"
)

// CmaEsCholB is optimize.CmaEsChol with xmin,xmax constraints
// only sendTask,ensureBounds are different
type CmaEsCholB struct {
	//optimize.CmaEsChol
	// InitStepSize sets the initial size of the covariance matrix adaptation.
	// If InitStepSize is 0, a default value of 0.5 is used. InitStepSize cannot
	// be negative, or CmaEsCholB will panic.
	InitStepSize float64
	// Population sets the population size for the algorithm. If Population is
	// 0, a default value of 4 + math.Floor(3*math.Log(float64(dim))) is used.
	// Population cannot be negative or CmaEsCholB will panic.
	Population int
	// InitCholesky specifies the Cholesky decomposition of the covariance
	// matrix for the initial sampling distribution. If InitCholesky is nil,
	// a default value of I is used. If it is non-nil, then it must have
	// InitCholesky.Size() be equal to the problem dimension.
	InitCholesky *mat.Cholesky
	// StopLogDet sets the threshold for stopping the optimization if the
	// distribution becomes too peaked. The log determinant is a measure of the
	// (log) "volume" of the normal distribution, and when it is too small
	// the samples are almost the same. If the log determinant of the covariance
	// matrix becomes less than StopLogDet, the optimization run is concluded.
	// If StopLogDet is 0, a default value of dim*log(1e-16) is used.
	// If StopLogDet is NaN, the stopping criterion is not used, though
	// this can cause numeric instabilities in the algorithm.
	StopLogDet float64
	// ForgetBest, when true, does not track the best overall function value found,
	// instead returning the new best sample in each iteration. If ForgetBest
	// is false, then the minimum value returned will be the lowest across all
	// iterations, regardless of when that sample was generated.
	ForgetBest bool
	// Src allows a random number generator to be supplied for generating samples.
	// If Src is nil the generator in golang.org/x/math/rand is used.
	Src rand.Source

	// Fixed algorithm parameters.
	dim                 int
	pop                 int
	weights             []float64
	muEff               float64
	cc, cs, c1, cmu, ds float64
	eChi                float64

	// Function data.
	xs *mat.Dense
	fs []float64

	// Adaptive algorithm parameters.
	invSigma float64 // inverse of the sigma parameter
	pc, ps   []float64
	mean     []float64
	chol     mat.Cholesky

	// Overall best.
	bestX, Xmin, Xmax []float64
	bestF             float64

	// Synchronization.
	sentIdx     int
	receivedIdx int
	operation   chan<- optimize.Task
	updateErr   error
}

var (
	_ optimize.Statuser = (*CmaEsCholB)(nil)
	_ optimize.Method   = (*CmaEsCholB)(nil)
)

// Needs ...
func (cma *CmaEsCholB) Needs() struct{ Gradient, Hessian bool } {
	return struct{ Gradient, Hessian bool }{false, false}
}

// Uses ...
func (cma *CmaEsCholB) Uses(has optimize.Available) (optimize.Available, error) {
	return optimize.Available{}, nil
}

func (cma *CmaEsCholB) methodConverged() optimize.Status {
	sd := cma.StopLogDet
	switch {
	case math.IsNaN(sd):
		return optimize.NotTerminated
	case sd == 0:
		sd = float64(cma.dim) * -36.8413614879 // ln(1e-16)
	}
	if cma.chol.LogDet() < sd {
		return optimize.MethodConverge
	}
	return optimize.NotTerminated
}

// Status returns the status of the method.
func (cma *CmaEsCholB) Status() (optimize.Status, error) {
	if cma.updateErr != nil {
		return optimize.Failure, cma.updateErr
	}
	return cma.methodConverged(), nil
}

// Init ...
func (cma *CmaEsCholB) Init(dim, tasks int) int {
	if dim <= 0 {
		panic(nonpositiveDimension)
	}
	if tasks < 0 {
		panic(negativeTasks)
	}

	// Set fixed algorithm parameters.
	// Parameter values are from https://arxiv.org/pdf/1604.00772.pdf .
	cma.dim = dim
	cma.pop = cma.Population
	n := float64(dim)
	if cma.pop == 0 {
		cma.pop = 4 + int(3*math.Log(n)) // Note the implicit floor.
	} else if cma.pop < 0 {
		panic("cma-es-chol: negative population size")
	}
	mu := cma.pop / 2
	cma.weights = resize(cma.weights, mu)
	for i := range cma.weights {
		v := math.Log(float64(mu)+0.5) - math.Log(float64(i)+1)
		cma.weights[i] = v
	}
	floats.Scale(1/floats.Sum(cma.weights), cma.weights)
	cma.muEff = 0
	for _, v := range cma.weights {
		cma.muEff += v * v
	}
	cma.muEff = 1 / cma.muEff

	cma.cc = (4 + cma.muEff/n) / (n + 4 + 2*cma.muEff/n)
	cma.cs = (cma.muEff + 2) / (n + cma.muEff + 5)
	cma.c1 = 2 / ((n+1.3)*(n+1.3) + cma.muEff)
	cma.cmu = math.Min(1-cma.c1, 2*(cma.muEff-2+1/cma.muEff)/((n+2)*(n+2)+cma.muEff))
	cma.ds = 1 + 2*math.Max(0, math.Sqrt((cma.muEff-1)/(n+1))-1) + cma.cs
	// E[chi] is taken from https://en.wikipedia.org/wiki/CMA-ES (there
	// listed as E[||N(0,1)||]).
	cma.eChi = math.Sqrt(n) * (1 - 1.0/(4*n) + 1/(21*n*n))

	// Allocate memory for function data.
	cma.xs = mat.NewDense(cma.pop, dim, nil)
	cma.fs = resize(cma.fs, cma.pop)

	// Allocate and initialize adaptive parameters.
	cma.invSigma = 1 / cma.InitStepSize
	if cma.InitStepSize == 0 {
		cma.invSigma = 10.0 / 3
	} else if cma.InitStepSize < 0 {
		panic("cma-es-chol: negative initial step size")
	}
	cma.pc = resize(cma.pc, dim)
	for i := range cma.pc {
		cma.pc[i] = 0
	}
	cma.ps = resize(cma.ps, dim)
	for i := range cma.ps {
		cma.ps[i] = 0
	}
	cma.mean = resize(cma.mean, dim) // mean location initialized at the start of Run

	if cma.InitCholesky != nil {
		if cma.InitCholesky.SymmetricDim() != dim {
			panic("cma-es-chol: incorrect InitCholesky size")
		}
		cma.chol.Clone(cma.InitCholesky)
	} else {
		// Set the initial Cholesky to I.
		b := mat.NewDiagDense(dim, nil)
		for i := 0; i < dim; i++ {
			b.SetDiag(i, 1)
		}
		var chol mat.Cholesky
		ok := chol.Factorize(b)
		if !ok {
			panic("cma-es-chol: bad cholesky. shouldn't happen")
		}
		cma.chol = chol
	}

	cma.bestX = resize(cma.bestX, dim)
	cma.bestF = math.Inf(1)

	cma.sentIdx = 0
	cma.receivedIdx = 0
	cma.operation = nil
	cma.updateErr = nil
	t := min(tasks, cma.pop)
	return t
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func (cma *CmaEsCholB) sendInitTasks(tasks []optimize.Task) {
	for i, task := range tasks {
		cma.sendTask(i, task)
	}
	cma.sentIdx = len(tasks)
}

func (cma *CmaEsCholB) ensureBounds(x []float64) {
	nBounded := 0
	for i := range x {
		if (i < len(cma.Xmin) && x[i] <= cma.Xmin[i]) || (i < len(cma.Xmax) && x[i] >= cma.Xmax[i]) {
			nBounded++
		}
	}
	for i := range x {
		if i < len(cma.Xmin) && x[i] < cma.Xmin[i] {
			if nBounded < len(x) {
				x[i] = cma.Xmin[i]
			} else {
				for x[i] < cma.Xmin[i] {
					x[i] = (x[i] + cma.mean[i]) / 2
				}
			}
		}
		if i < len(cma.Xmax) && x[i] > cma.Xmax[i] {
			if nBounded < len(x) {
				x[i] = cma.Xmax[i]
			} else {
				for x[i] > cma.Xmax[i] {
					x[i] = (x[i] + cma.mean[i]) / 2
				}
			}
		}
	}
}

// sendTask generates a sample and sends the task. It does not update the cma index.
// this method differs of original cmaes in using ensureBounds
func (cma *CmaEsCholB) sendTask(idx int, task optimize.Task) {
	task.ID = idx
	task.Op = optimize.FuncEvaluation
	distmv.NormalRand(cma.xs.RawRowView(idx), cma.mean, &cma.chol, cma.Src)
	cma.ensureBounds(cma.xs.RawRowView(idx))
	copy(task.X, cma.xs.RawRowView(idx))
	cma.operation <- task
}

// bestIdx returns the best index in the functions. Returns -1 if all values
// are NaN.
func (cma *CmaEsCholB) bestIdx() int {
	best := -1
	bestVal := math.Inf(1)
	for i, v := range cma.fs {
		if math.IsNaN(v) {
			continue
		}
		// Use equality in case somewhere evaluates to +inf.
		if v <= bestVal {
			best = i
			bestVal = v
		}
	}
	return best
}

// findBestAndUpdateTask finds the best task in the current list, updates the
// new best overall, and then stores the best location into task.
func (cma *CmaEsCholB) findBestAndUpdateTask(task optimize.Task) optimize.Task {
	// Find and update the best location.
	// Don't use floats because there may be NaN values.
	best := cma.bestIdx()
	bestF := math.NaN()
	bestX := cma.xs.RawRowView(0)
	if best != -1 {
		bestF = cma.fs[best]
		bestX = cma.xs.RawRowView(best)
	}
	if cma.ForgetBest {
		task.F = bestF
		copy(task.X, bestX)
	} else {
		if bestF < cma.bestF {
			cma.bestF = bestF
			copy(cma.bestX, bestX)
		}
		task.F = cma.bestF
		copy(task.X, cma.bestX)
	}
	return task
}

// Run ...
func (cma *CmaEsCholB) Run(operations chan<- optimize.Task, results <-chan optimize.Task, tasks []optimize.Task) {
	copy(cma.mean, tasks[0].X)
	cma.operation = operations
	// Send the initial tasks. We know there are at most as many tasks as elements
	// of the population.
	cma.sendInitTasks(tasks)

Loop:
	for {
		result := <-results
		switch result.Op {
		default:
			panic("unknown operation")
		case optimize.PostIteration:
			break Loop
		case optimize.MajorIteration:
			// The last thing we did was update all of the tasks and send the
			// major iteration. Now we can send a group of tasks again.
			cma.sendInitTasks(tasks)
		case optimize.FuncEvaluation:
			cma.receivedIdx++
			cma.fs[result.ID] = result.F
			switch {
			case cma.sentIdx < cma.pop:
				// There are still tasks to evaluate. Send the next.
				cma.sendTask(cma.sentIdx, result)
				cma.sentIdx++
			case cma.receivedIdx < cma.pop:
				// All the tasks have been sent, but not all of them have been received.
				// Need to wait until all are back.
				continue Loop
			default:
				// All of the evaluations have been received.
				if cma.receivedIdx != cma.pop {
					panic("bad logic")
				}
				cma.receivedIdx = 0
				cma.sentIdx = 0

				task := cma.findBestAndUpdateTask(result)
				// Update the parameters and send a MajorIteration or a convergence.
				err := cma.update()
				// Kill the existing data.
				for i := range cma.fs {
					cma.fs[i] = math.NaN()
					cma.xs.Set(i, 0, math.NaN())
				}
				switch {
				case err != nil:
					cma.updateErr = err
					task.Op = optimize.MethodDone
				case cma.methodConverged() != optimize.NotTerminated:
					task.Op = optimize.MethodDone
				default:
					task.Op = optimize.MajorIteration
					task.ID = -1
				}
				operations <- task
			}
		}
	}

	// Been told to stop. Clean up.
	// Need to see best of our evaluated tasks so far. Should instead just
	// collect, then see.
	for task := range results {
		switch task.Op {
		case optimize.MajorIteration:
		case optimize.FuncEvaluation:
			cma.fs[task.ID] = task.F
		default:
			panic("unknown operation")
		}
	}
	// Send the new best value if the evaluation is better than any we've
	// found so far. Keep this separate from findBestAndUpdateTask so that
	// we only send an iteration if we find a better location.
	if !cma.ForgetBest {
		best := cma.bestIdx()
		if best != -1 && cma.fs[best] < cma.bestF {
			task := tasks[0]
			task.F = cma.fs[best]
			copy(task.X, cma.xs.RawRowView(best))
			task.Op = optimize.MajorIteration
			task.ID = -1
			operations <- task
		}
	}
	close(operations)
}

// update computes the new parameters (mean, cholesky, etc.). Does not update
// any of the synchronization parameters (taskIdx).
func (cma *CmaEsCholB) update() error {
	// Sort the function values to find the elite samples.
	ftmp := make([]float64, cma.pop)
	copy(ftmp, cma.fs)
	indexes := make([]int, cma.pop)
	for i := range indexes {
		indexes[i] = i
	}
	sort.Sort(bestSorter{F: ftmp, Idx: indexes})

	meanOld := make([]float64, len(cma.mean))
	copy(meanOld, cma.mean)

	// m_{t+1} = \sum_{i=1}^mu w_i x_i
	for i := range cma.mean {
		cma.mean[i] = 0
	}
	for i, w := range cma.weights {
		idx := indexes[i] // index of teh 1337 sample.
		floats.AddScaled(cma.mean, w, cma.xs.RawRowView(idx))
	}
	cma.ensureBounds(cma.mean)
	meanDiff := make([]float64, len(cma.mean))
	floats.SubTo(meanDiff, cma.mean, meanOld)

	// p_{c,t+1} = (1-c_c) p_{c,t} + \sqrt(c_c*(2-c_c)*mueff) (m_{t+1}-m_t)/sigma_t
	floats.Scale(1-cma.cc, cma.pc)
	scaleC := math.Sqrt(cma.cc*(2-cma.cc)*cma.muEff) * cma.invSigma
	floats.AddScaled(cma.pc, scaleC, meanDiff)

	// p_{sigma, t+1} = (1-c_sigma) p_{sigma,t} + \sqrt(c_s*(2-c_s)*mueff) A_t^-1 (m_{t+1}-m_t)/sigma_t
	floats.Scale(1-cma.cs, cma.ps)
	// First compute A_t^-1 (m_{t+1}-m_t), then add the scaled vector.
	tmp := make([]float64, cma.dim)
	tmpVec := mat.NewVecDense(cma.dim, tmp)
	diffVec := mat.NewVecDense(cma.dim, meanDiff)
	err := tmpVec.SolveVec(cma.chol.RawU().T(), diffVec)
	if err != nil {
		return err
	}
	scaleS := math.Sqrt(cma.cs*(2-cma.cs)*cma.muEff) * cma.invSigma
	floats.AddScaled(cma.ps, scaleS, tmp)

	// Compute the update to A.
	scaleChol := 1 - cma.c1 - cma.cmu
	if scaleChol == 0 {
		scaleChol = math.SmallestNonzeroFloat64 // enough to kill the old data, but still non-zero.
	}
	cma.chol.Scale(scaleChol, &cma.chol)
	cma.chol.SymRankOne(&cma.chol, cma.c1, mat.NewVecDense(cma.dim, cma.pc))
	for i, w := range cma.weights {
		idx := indexes[i]
		floats.SubTo(tmp, cma.xs.RawRowView(idx), meanOld)
		cma.chol.SymRankOne(&cma.chol, cma.cmu*w*cma.invSigma, tmpVec)
	}

	// sigma_{t+1} = sigma_t exp(c_sigma/d_sigma * norm(p_{sigma,t+1}/ E[chi] -1)
	normPs := floats.Norm(cma.ps, 2)
	cma.invSigma /= math.Exp(cma.cs / cma.ds * (normPs/cma.eChi - 1))
	return nil
}

type bestSorter struct {
	F   []float64
	Idx []int
}

func (b bestSorter) Len() int {
	return len(b.F)
}
func (b bestSorter) Less(i, j int) bool {
	return b.F[i] < b.F[j]
}
func (b bestSorter) Swap(i, j int) {
	b.F[i], b.F[j] = b.F[j], b.F[i]
	b.Idx[i], b.Idx[j] = b.Idx[j], b.Idx[i]
}
