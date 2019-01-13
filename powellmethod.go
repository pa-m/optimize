package optimize

// Powell optimize.Method
/*
var (
	_ optimize.Statuser = (*Powell)(nil)
	_ optimize.Method   = (*Powell)(nil)
)

// Powell optimize.Method
type Powell struct {
	status optimize.Status
	err    error
}

// Needs for Powell
func (m *Powell) Needs() struct{ Gradient, Hessian bool } {
	return struct{ Gradient, Hessian bool }{false, false}
}

// Status for Powell
func (m *Powell) Status() (optimize.Status, error) {
	var s optimize.Status
	return s, nil
}

// Init for Powell
func (m *Powell) Init(dim, tasks int) int {
	m.status = optimize.NotTerminated
	m.err = nil
	return 0
}

// Run for Powell
func (m *Powell) Run(operations chan<- optimize.Task, results <-chan optimize.Task, tasks []optimize.Task) {
	for task := range results {
		switch task.Op {
		}
	}
	close(operations)
}
*/
