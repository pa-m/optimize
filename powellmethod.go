// Copyright ©2016 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package optimize

import (
	"math"

	"gonum.org/v1/gonum/optimize"
)

// Powell is a global optimizer that evaluates the function at random
// locations. Not a good optimizer, but useful for comparison and debugging.
type Powell struct {
	PM       *PowellMinimizer
	settings optimize.Settings
	status   optimize.Status
	err      error
	bestF    float64
	bestX    []float64
}

// Needs for Powell to implement gonum optimize.Needser
func (g *Powell) Needs() struct{ Gradient, Hessian bool } {
	return struct{ Gradient, Hessian bool }{false, false}
}

// Init for Powell to implement gonum optimize.Method
func (g *Powell) Init(dim, tasks int) int {
	if dim <= 0 {
		panic(nonpositiveDimension)
	}
	if tasks < 0 {
		panic(negativeTasks)
	}
	g.bestF = math.Inf(1)
	g.bestX = resize(g.bestX, dim)
	return 1
}

func (g *Powell) updateMajor(operation chan<- optimize.Task, task optimize.Task) {
	// Update the best value seen so far, and send a MajorIteration.
	if task.F < g.bestF {
		g.bestF = task.F
		copy(g.bestX, task.X)
	}
	task.Op = optimize.MajorIteration
	operation <- task
}

// Run for Powell to implement gonum optimize.Method
func (g *Powell) Run(operation chan<- optimize.Task, result <-chan optimize.Task, tasks []optimize.Task) {
	var stop bool
	fnMaxIter := func(int) bool { return stop }
	fnMaxFev := func(int) bool { return stop }

	if g.PM == nil {
		g.PM = NewPowellMinimizer()
	}
	pm := g.PM

	result1 := make(chan optimize.Task)
	// Send initial tasks to evaluate

	dup := func(x []float64) []float64 {
		r := make([]float64, len(x))
		copy(r, x)
		return r
	}
	InitX := tasks[0].Location.X
	go func(id int) {
		_, warnflag := minimizePowell(func(x []float64) (y float64) {
			y = math.NaN()
			defer func() {
				if r := recover(); r == "send on closed channel" {
					return
				}
			}()
			operation <- optimize.Task{ID: id, Op: optimize.FuncEvaluation, Location: &optimize.Location{X: dup(x)}}
			task := <-result1
			if task.Location != nil {
				y = task.Location.F
			}
			return
		}, InitX, nil, pm.Xtol, pm.Ftol, fnMaxIter, fnMaxFev, pm.Logger)
		switch warnflag {
		case 1:
			g.status = optimize.FunctionEvaluationLimit
		case 2:
			g.status = optimize.IterationLimit
		default:
			g.status = optimize.MethodConverge
		}

		defer func() {
			if r := recover(); r == "send on closed channel" {
				return
			}
		}()
		operation <- optimize.Task{ID: id, Op: optimize.MethodDone}

	}(0)

	// Read from the channel until PostIteration is sent.
Loop:
	for {
		task := <-result
		switch task.Op {
		default:
			panic("unknown operation")
		case optimize.NoOperation, optimize.PostIteration:
			close(result1)
			break Loop
		case optimize.MajorIteration:

		case optimize.FuncEvaluation:
			result1 <- task
			g.updateMajor(operation, task)
		}
	}

	// PostIteration was sent. Update the best new values.
	for task := range result {
		switch task.Op {
		default:
			panic("unknown operation")
		case optimize.MajorIteration:
		case optimize.FuncEvaluation:
			g.updateMajor(operation, task)
		case optimize.NoOperation:
		}
	}
	stop = true
	close(operation)
}

// Status ...
func (g *Powell) Status() (optimize.Status, error) {
	return g.status, g.err
}
