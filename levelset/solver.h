/*
 *  solver.h
 *  smoke
 *
 */

#include "common.h"

namespace solver {
	// Solve Ax = b
	FLOAT solve( FLOAT **A, FLOAT **x, FLOAT **b, int n, char subcell, char solver );
}
