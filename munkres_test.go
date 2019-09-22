package munkres_test

import (
	"fmt"
	"math"
	"reflect"
	"testing"

	"github.com/charles-haynes/munkres"
)

func computeCost(matrix [][]float64, match []int) (float64, error) {
	result := 0.0
	visited := map[int]int{}
	for i := 0; i < len(matrix); i++ {
		if match[i] == -1 {
			continue
		}
		if _, ok := visited[match[i]]; ok {
			return 0.0, fmt.Errorf(
				"workers %d and %d have the same job",
				match[i], i)
		}
		visited[match[i]] = i
		result += matrix[i][match[i]]
	}
	return result, nil
}

type test struct {
	name       string
	costMatrix [][]float64
	err        error
	res        []int
	cost       float64
}

var tests = []test{
	test{
		"test1",
		[][]float64{
			[]float64{4.0, 1.5, 4.0},
			[]float64{4.0, 4.5, 6.0},
			[]float64{3.0, 2.25, 3.0},
		},
		nil,
		[]int{1, 0, 2},
		1.5 + 4.0 + 3.0,
	},
	test{
		"test2",
		[][]float64{
			[]float64{1.0, 1.0, 0.8},
			[]float64{0.9, 0.8, 0.1},
			[]float64{0.9, 0.7, 0.4},
		},
		nil,
		[]int{0, 2, 1},
		1.0 + 0.1 + 0.7,
	},
	test{
		"test3",
		[][]float64{
			[]float64{6.0, 0.0, 7.0, 5.0},
			[]float64{2.0, 6.0, 2.0, 6.0},
			[]float64{2.0, 7.0, 2.0, 1.0},
			[]float64{9.0, 4.0, 7.0, 1.0},
		},
		nil,
		[]int{1, 0, 2, 3},
		0.0 + 2.0 + 2.0 + 1.0,
	},
	test{
		"not rectangular",
		[][]float64{
			[]float64{1.0, 2.0},
			[]float64{3.0},
		},
		munkres.ErrorIrregularCostMatrix,
		nil,
		0.0,
	},
	test{
		"infinite cost",
		[][]float64{
			[]float64{1.0, 2.0},
			[]float64{3.0, math.Inf(1)},
		},
		munkres.ErrorInfiniteCost,
		nil,
		0.0,
	},
	test{
		"NaN cost",
		[][]float64{
			[]float64{1.0, 2.0},
			[]float64{3.0, math.NaN()},
		},
		munkres.ErrorNaNCost,
		nil,
		0.0,
	},
	test{
		"nil matrix",
		nil,
		nil,
		nil,
		0.0,
	},
	test{
		"empty matrix",
		[][]float64{{}},
		nil,
		[]int{-1},
		0.0,
	},
	test{
		"unassigned job",
		[][]float64{
			[]float64{6.0, 0.0, 7.0, 5.0, 2.0},
			[]float64{2.0, 6.0, 2.0, 6.0, 7.0},
			[]float64{2.0, 7.0, 2.0, 1.0, 1.0},
			[]float64{9.0, 4.0, 7.0, 1.0, 0.0},
		},
		nil,
		[]int{1, 0, 3, 4},
		0.0 + 2.0 + 1.0 + 0.0,
	},
	test{
		"unassigned worker",
		[][]float64{
			[]float64{6.0, 0.0, 7.0, 5.0},
			[]float64{2.0, 6.0, 2.0, 6.0},
			[]float64{2.0, 7.0, 2.0, 1.0},
			[]float64{9.0, 4.0, 7.0, 1.0},
			[]float64{0.0, 0.0, 0.0, 0.0},
		},
		nil,
		[]int{1, -1, 2, 3, 0},
		0.0 + 2.0 + 1.0 + 0.0,
	},
}

// create a complete nxn matrix test case where each cell i,j cost is (i+1)*(j+1)
// the solution is the diagonal from (0,n-1)..(n-1,0)
func CreateTest(n int) test {
	test := test{
		name: fmt.Sprintf("%dx%d matrix", n, n),
	}
	test.costMatrix = make([][]float64, n)
	for i := 0; i < n; i++ {
		test.costMatrix[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			test.costMatrix[i][j] = float64((i + 1) * (j + 1))
		}
	}
	test.res = make([]int, n)
	for i := 0; i < n; i++ {
		test.res[i] = n - i - 1
		test.cost += float64((i + 1) * (n - i))
	}
	return test
}

func TestAbs(t *testing.T) {
	tests := append(tests, CreateTest(100))
	for _, d := range tests {
		h, err := munkres.NewHungarianAlgorithm(d.costMatrix)
		if err != d.err {
			t.Errorf("%s: want err = %s got %s",
				d.name, d.err, err)
		}
		if d.err != nil {
			continue
		}
		res := h.Execute()
		if !reflect.DeepEqual(res, d.res) {
			t.Errorf("%s: want res = %v got %v",
				d.name, d.res, res)
		}
		cost, err := computeCost(d.costMatrix, res)
		if err != nil {
			t.Errorf("%s: %s", d.name, err)
		}
		if math.Abs(cost-d.cost) > 0.0000001 {
			t.Errorf("%s: want cost = %f got %f",
				d.name, d.cost, cost)
		}
	}
}
