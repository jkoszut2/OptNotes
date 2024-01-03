// Math Functions

#include "prob_workspace.h"

double coneprod1(double u, double v, int n) {
	// Cone product of vectors u, v in an order-1 cone
	// n is the number of elements in u or v
	double uv = 0;
	for (int i=0; i<n; i++) {
	    uv += cone_orders[i] * cone_exprs[i];
	}
	
	return uv;
}

int coneprod2(double *u, double *v, int n, double *uv) {
	// Cone product of vectors u, v with n elements and
	// belonging to a cone with order >1
	
	for (int i=0; i<n; i++) {
	    uv[0] += u[i] * v[i];
	}
	for (int i=1; i<n; i++) {
	    uv[i] = u[0]*v[i] + v[0]*u[i];
	}

	return 0;
}