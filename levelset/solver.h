/*
 *  solver.h
 *  smoke
 *
 */

#include "common.h"
#include <CL/cl.h>

namespace solver {
	// Solve Ax = b
	FLOAT solve( FLOAT ***A, FLOAT ***x, FLOAT ***b, char subcell, char solver );
	void setCL(int n);
}
