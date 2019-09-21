An implementation of the _Hungarian algorithm_ for solving the assignment
problem. An instance of the assignment problem consists of a number of
workers along with a number of jobs and a cost matrix which gives the cost of
assigning the i'th worker to the j'th job at position (i, j). The goal is to
find an assignment of workers to jobs so that no job is assigned more than
one worker and so that no worker is assigned to more than one job in such a
manner so as to minimize the total cost of completing the jobs.

An assignment for a cost matrix that has more workers than jobs will
necessarily include unassigned workers, indicated by an assignment value of
-1; in no other circumstance will there be unassigned workers. Similarly, an
assignment for a cost matrix that has more jobs than workers will necessarily
include unassigned jobs; in no other circumstance will there be unassigned
jobs. For completeness, an assignment for a square cost matrix will give
exactly one unique worker to each job.

This version of the Hungarian algorithm runs in time O(n^3), where n is the
maximum among the number of workers and the number of jobs.

This is a Go implementation of the Java version found at
https://github.com/KevinStern/software-and-algorithms/
