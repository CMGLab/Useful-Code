void mapToStandardInterval(double & x_std, double x, double old_min,
													 double old_max, double new_min, double new_max);
void quadratureData(double * quad_pts, double * quad_wts, int num_quad_pts,
										double interval_start, double interval_end);

void quadratureData(double * quad_pts, double * quad_wts, int num_quad_pts,
										double interval_start, double interval_end)
/*
Generates the Gaussian quadrature scheme based on the Legendre nodes containing
N points and weights which is accurate up to degree p = 2*N-1. The interval of
integration is [a,b]. If a and be are not specified by the user, the
default interval is [0,1].

   Inputs: num_quad_pts - Number of quadrature points and weights desired
           interval_start - left endpoint of integration interval
           interval_end - right endpoint of integration interval

   Outputs: quad_pts - list of N quadrature nodes
            quad_wts - list of N quadrature weights

Based on the algorithm introduced in:
   Golub, Gene H., and John H. Welsch. "Calculation of Gauss quadrature rules."
   Mathematics of computation 23.106 (1969): 221-230.

 Written by: Joseph Benzaken
 Last Updated: 5/4/16

 CMGLab
*/
{
	double * beta = new double[num_quad_pts];
		for (int i = 0; i < num_quad_pts; i++)
			beta[i] = 0.0;

	//Legendre 3-term recursion relationship
	for (int i = 0; i < num_quad_pts; i++) {
		int j = i + 1;
		beta[i] = (double) (j-1)*(j-1) / (double) (4*(j-1)*(j-1)-1);
	}

	double * J = new double[num_quad_pts*num_quad_pts];
	for (int i = 0; i < num_quad_pts*num_quad_pts; i++)
		J[i] = 0.0;
	//assemble Jacobi matrix
	for (int i = 1; i < num_quad_pts; i++) {
		J[(i-1)*num_quad_pts + i] = sqrt(beta[i]);
		J[i*num_quad_pts + i - 1] = sqrt(beta[i]);
	}
	//compute eigenvalues and eigenvectors of Jacobi matrix
	LAPACKE_dsyev(CblasRowMajor, 'V', 'U', num_quad_pts, J, num_quad_pts, quad_pts);

	for (int i = 0; i < num_quad_pts; i++) {
		mapToStandardInterval(quad_pts[i], quad_pts[i], -1.0, 1.0, interval_start,
													interval_end);
		quad_wts[i] = J[i]*J[i]*(interval_end-interval_start);
	}

	delete [] beta;
	delete [] J;
}

void mapToStandardInterval(double & x_std, double x, double old_min,
													 double old_max, double new_min, double new_max)
{
	x_std = (new_max - new_min) / (old_max - old_min) * (x - old_min) + new_min;
}
