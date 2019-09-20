package munkres

import (
	"fmt"
	"math"
)

/*
 * An implementation of the Hungarian algorithm for solving the assignment
 * problem. An instance of the assignment problem consists of a number of
 * workers along with a number of jobs and a cost matrix which gives the cost of
 * assigning the i'th worker to the j'th job at position (i, j). The goal is to
 * find an assignment of workers to jobs so that no job is assigned more than
 * one worker and so that no worker is assigned to more than one job in such a
 * manner so as to minimize the total cost of completing the jobs.
  *
 * An assignment for a cost matrix that has more workers than jobs will
 * necessarily include unassigned workers, indicated by an assignment value of
 * -1; in no other circumstance will there be unassigned workers. Similarly, an
 * assignment for a cost matrix that has more jobs than workers will necessarily
 * include unassigned jobs; in no other circumstance will there be unassigned
 * jobs. For completeness, an assignment for a square cost matrix will give
 * exactly one unique worker to each job.
 *
 * This version of the Hungarian algorithm runs in time O(n^3), where n is the
 * maximum among the number of workers and the number of jobs.
 *
 * ported from the Java version by Kevin L. Stern
*/
type HungarianAlgorithm struct {
	costMatrix                         [][]float64
	rows, cols, dim                    int
	labelByWorker, labelByJob          []float64
	minSlackWorkerByJob                []int
	minSlackValueByJob                 []float64
	matchJobByWorker, matchWorkerByJob []int
	parentWorkerByCommittedJob         []int
	committedWorkers                   []bool
}

// Construct an instance of the algorithm.
//
// @param costMatrix
//          the cost matrix, where matrix[i][j] holds the cost of assigning
//          worker i to job j, for all i, j. The cost matrix must not be
//          irregular in the sense that all rows must be the same length; in
//          addition, all entries must be non-infinite numbers.
func NewHungarianAlgorithm(costMatrix [][]float64) (HungarianAlgorithm, error) {
	dim := len(costMatrix)
	if len(costMatrix[0]) > dim {
		dim = len(costMatrix[0])
	}
	this := HungarianAlgorithm{
		costMatrix:                 make([][]float64, dim),
		rows:                       len(costMatrix),
		cols:                       len(costMatrix[0]),
		dim:                        dim,
		labelByWorker:              make([]float64, dim),
		labelByJob:                 make([]float64, dim),
		minSlackWorkerByJob:        make([]int, dim),
		minSlackValueByJob:         make([]float64, dim),
		committedWorkers:           make([]bool, dim),
		parentWorkerByCommittedJob: make([]int, dim),
		matchJobByWorker:           make([]int, dim),
		matchWorkerByJob:           make([]int, dim),
	}
	for w := 0; w < dim; w++ {
		this.costMatrix[w] = make([]float64, dim)
		if w > len(costMatrix) {
			continue
		}
		if len(costMatrix[w]) != this.cols {
			return this, fmt.Errorf("Irregular cost matrix")
		}
		for j := range costMatrix[w] {
			if math.IsInf(costMatrix[w][j], 0) {
				return this, fmt.Errorf("Infinite cost")
			}
			if math.IsNaN(costMatrix[w][j]) {
				return this, fmt.Errorf("NaN cost")
			}
		}
		copy(this.costMatrix[w], costMatrix[w])
	}
	for i := 0; i < dim; i++ {
		this.matchJobByWorker[i] = -1
		this.matchWorkerByJob[i] = -1
	}
	return this, nil
}

// Compute an initial feasible solution by assigning zero labels to the
// workers and by assigning to each job a label equal to the minimum cost
// among its incident edges.
func (h *HungarianAlgorithm) computeInitialFeasibleSolution() {
	for j := range h.labelByJob {
		h.labelByJob[j] = math.Inf(1)
	}
	for w := 0; w < h.dim; w++ {
		for j := 0; j < h.dim; j++ {
			if h.costMatrix[w][j] < h.labelByJob[j] {
				h.labelByJob[j] = h.costMatrix[w][j]
			}
		}
	}
}

// Execute the algorithm.
//
// @return the minimum cost matching of workers to jobs based upon the
//         provided cost matrix. A matching value of -1 indicates that the
//         corresponding worker is unassigned.
func (h *HungarianAlgorithm) Execute() []int {
	// Heuristics to improve performance: Reduce rows and columns by their
	// smallest element, compute an initial non-zero dual feasible solution and
	// create a greedy matching from workers to jobs of the cost matrix.
	h.reduce()
	h.computeInitialFeasibleSolution()
	h.greedyMatch()

	for w := h.fetchUnmatchedWorker(); w < h.dim; w = h.fetchUnmatchedWorker() {
		h.initializePhase(w)
		h.executePhase()
	}
	result := h.matchJobByWorker[:h.rows]
	for w := range result {
		if result[w] >= h.cols {
			result[w] = -1
		}
	}
	return result
}

// Execute a single phase of the algorithm. A phase of the Hungarian algorithm
// consists of building a set of committed workers and a set of committed jobs
// from a root unmatched worker by following alternating unmatched/matched
// zero-slack edges. If an unmatched job is encountered, then an augmenting
// path has been found and the matching is grown. If the connected zero-slack
// edges have been exhausted, the labels of committed workers are increased by
// the minimum slack among committed workers and non-committed jobs to create
// more zero-slack edges (the labels of committed jobs are simultaneously
// decreased by the same amount in order to maintain a feasible labeling).
// <p>
//
// The runtime of a single phase of the algorithm is O(n^2), where n is the
// dimension of the internal square cost matrix, since each edge is visited at
// most once and since increasing the labeling is accomplished in time O(n) by
// maintaining the minimum slack values among non-committed jobs. When a phase
// completes, the matching will have increased in size.
func (h *HungarianAlgorithm) executePhase() {
	for {
		minSlackWorker := -1
		minSlackJob := -1
		minSlackValue := math.Inf(1)
		for j := 0; j < h.dim; j++ {
			if h.parentWorkerByCommittedJob[j] == -1 {
				if h.minSlackValueByJob[j] < minSlackValue {
					minSlackValue = h.minSlackValueByJob[j]
					minSlackWorker = h.minSlackWorkerByJob[j]
					minSlackJob = j
				}
			}
		}
		if minSlackValue > 0 {
			h.updateLabeling(minSlackValue)
		}
		h.parentWorkerByCommittedJob[minSlackJob] = minSlackWorker
		if h.matchWorkerByJob[minSlackJob] == -1 {
			// An augmenting path has been found.
			committedJob := minSlackJob
			parentWorker := h.parentWorkerByCommittedJob[committedJob]
			for {
				temp := h.matchJobByWorker[parentWorker]
				h.match(parentWorker, committedJob)
				committedJob = temp
				if committedJob == -1 {
					break
				}
				parentWorker = h.parentWorkerByCommittedJob[committedJob]
			}
			return
		} else {
			// Update slack values since we increased the size of the committed
			// workers set.
			worker := h.matchWorkerByJob[minSlackJob]
			h.committedWorkers[worker] = true
			for j := 0; j < h.dim; j++ {
				if h.parentWorkerByCommittedJob[j] == -1 {
					slack := h.costMatrix[worker][j] -
						h.labelByWorker[worker] -
						h.labelByJob[j]
					if h.minSlackValueByJob[j] > slack {
						h.minSlackValueByJob[j] = slack
						h.minSlackWorkerByJob[j] = worker
					}
				}
			}
		}
	}
}

// @return the first unmatched worker or {@link #dim} if none.
func (h *HungarianAlgorithm) fetchUnmatchedWorker() int {
	for w, v := range h.matchJobByWorker {
		if v == -1 {
			return w
		}
	}
	return h.dim
}

// Find a valid matching by greedily selecting among zero-cost matchings. This
// is a heuristic to jump-start the augmentation algorithm.
func (h *HungarianAlgorithm) greedyMatch() {
	for w := 0; w < h.dim; w++ {
		for j := 0; j < h.dim; j++ {
			if h.matchJobByWorker[w] == -1 &&
				h.matchWorkerByJob[j] == -1 &&
				h.costMatrix[w][j]-h.labelByWorker[w]-
					h.labelByJob[j] == 0 {
				h.match(w, j)
			}
		}
	}
}

// Initialize the next phase of the algorithm by clearing the committed
// workers and jobs sets and by initializing the slack arrays to the values
// corresponding to the specified root worker.
//
// @param w
//          the worker at which to root the next phase.
func (h *HungarianAlgorithm) initializePhase(w int) {
	for i := range h.committedWorkers {
		h.committedWorkers[i] = false
	}
	for i := range h.parentWorkerByCommittedJob {
		h.parentWorkerByCommittedJob[i] = -1
	}
	h.committedWorkers[w] = true
	for j := 0; j < h.dim; j++ {
		h.minSlackValueByJob[j] = h.costMatrix[w][j] -
			h.labelByWorker[w] -
			h.labelByJob[j]
		h.minSlackWorkerByJob[j] = w
	}
}

// Helper method to record a matching between worker w and job j.
func (h *HungarianAlgorithm) match(w, j int) {
	h.matchJobByWorker[w] = j
	h.matchWorkerByJob[j] = w
}

// Reduce the cost matrix by subtracting the smallest element of each row from
// all elements of the row as well as the smallest element of each column from
// all elements of the column. Note that an optimal assignment for a reduced
// cost matrix is optimal for the original cost matrix.
func (h *HungarianAlgorithm) reduce() {
	for w := 0; w < h.dim; w++ {
		min := math.Inf(1)
		for j := 0; j < h.dim; j++ {
			if h.costMatrix[w][j] < min {
				min = h.costMatrix[w][j]
			}
		}
		for j := 0; j < h.dim; j++ {
			h.costMatrix[w][j] -= min
		}
	}
	min := make([]float64, h.dim)
	for j := 0; j < h.dim; j++ {
		min[j] = math.Inf(1)
	}
	for w := 0; w < h.dim; w++ {
		for j := 0; j < h.dim; j++ {
			if h.costMatrix[w][j] < min[j] {
				min[j] = h.costMatrix[w][j]
			}
		}
	}
	for w := 0; w < h.dim; w++ {
		for j := 0; j < h.dim; j++ {
			h.costMatrix[w][j] -= min[j]
		}
	}
}

// Update labels with the specified slack by adding the slack value for
// committed workers and by subtracting the slack value for committed jobs. In
// addition, update the minimum slack values appropriately.
func (h *HungarianAlgorithm) updateLabeling(slack float64) {
	for w := 0; w < h.dim; w++ {
		if h.committedWorkers[w] {
			h.labelByWorker[w] += slack
		}
	}
	for j := 0; j < h.dim; j++ {
		if h.parentWorkerByCommittedJob[j] != -1 {
			h.labelByJob[j] -= slack
		} else {
			h.minSlackValueByJob[j] -= slack
		}
	}
}

/* Example
func main() {
	k := 100
	fmt.Printf("Starting k = %d\n", k)
	start := time.Now()
	c := make([][]float64, k)
	for i := 0; i < k; i++ {
		c[i] = make([]float64, k)
		for j := 0; j < k; j++ {
			c[i][j] = rand.Float64()
		}
	}
	h, err := NewHungarianAlgorithm(c)
	if err != nil {
		panic(err)
	}
	r := h.Execute()
	fmt.Printf("Took: %s\nResult: %v\n", time.Since(start), r)
}
*/
