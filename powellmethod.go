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
	bestF float64
	bestX []float64
}

func (g *Powell) Needs() struct{ Gradient, Hessian bool } {
	return struct{ Gradient, Hessian bool }{false, false}
}

func (g *Powell) Init(dim, tasks int) int {
	if dim <= 0 {
		panic(nonpositiveDimension)
	}
	if tasks < 0 {
		panic(negativeTasks)
	}
	g.bestF = math.Inf(1)
	g.bestX = resize(g.bestX, dim)
	return tasks
}

func (g *Powell) sendNewLoc(operation chan<- optimize.Task, task optimize.Task) {
	// TODO set task X
	task.Op = optimize.FuncEvaluation
	operation <- task
}

func (g *Powell) updateMajor(operation chan<- optimize.Task, task optimize.Task) {
	// Update the best value seen so far, and send a MajorIteration.
	if task.F < g.bestF {
		g.bestF = task.F
		copy(g.bestX, task.X)
	} else {
		task.F = g.bestF
		copy(task.X, g.bestX)
	}
	task.Op = optimize.MajorIteration
	operation <- task
}

func (g *Powell) Run(operation chan<- optimize.Task, result <-chan optimize.Task, tasks []optimize.Task) {
	// Send initial tasks to evaluate
	for _, task := range tasks {
		g.sendNewLoc(operation, task)
	}

	// Read from the channel until PostIteration is sent.
Loop:
	for {
		task := <-result
		switch task.Op {
		default:
			panic("unknown operation")
		case optimize.PostIteration:
			break Loop
		case optimize.MajorIteration:
			g.sendNewLoc(operation, task)
		case optimize.FuncEvaluation:
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
		}
	}
	close(operation)
}
