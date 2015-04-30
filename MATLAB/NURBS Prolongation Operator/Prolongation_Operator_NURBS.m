function [P,new_w,new_n_1,new_Xi_1,new_n_2,new_Xi_2] = Prolongation_Operator_NURBS(p_1,n_1,Xi_1,add_Xi_1,p_2,n_2,Xi_2,add_Xi_2,w)

% Prolongation_Operator_NURBS
% Written by Christopher Coley
% Last modified 30 Apr 15
%
% Creates the prolongation operator for 2D with NURBS basis
%
% Output:
%   P - 2D prolongation operator
%   new_w - NURBS weights for the refined NURBS basis
%   new_n_1 - number of basis functions for the refined B-spline space in
%   first dimension
%   new_Xi_1 - univariate knot vector for the refined B-spline space in
%   first
%   dimension
%   new_n_2 - number of basis functions for the refined B-spline space in
%   second dimension
%   new_Xi_2 - univariate knot vector for the refined B-spline space in second
%   dimension
%
%
% Input:
%   p_1 - polynomial degree of the original and refined NURBS curves in
%   first dimension;  scalar
%   n_1 - number of basis functions for the original NURBS space in first
%   dimension; scalar
%   Xi_1 - univariate knot vector for the original NURBS space in first
%   dimension; row vector
%   add_Xi_1 - knots to be added in first dimension; row vector
%   p_2 - polynomial degree of the original and refined NURBS curves in
%   second dimension;  scalar
%   n_2 - number of basis functions for the original NURBS space in second
%   dimension; scalar
%   Xi_2 - univariate knot vector for the original NURBS space in second
%   dimension; row vector
%   add_Xi_2 - knots to be added in second dimension; row vector
%   w - NURBS weights for the unrefined space; 1 x (n_1 x n_2) vector;

Wc = diag(w);
[T,new_n_1,new_Xi_1,new_n_2,new_Xi_2] = Prolongation_Operator_2D(p_1,n_1,Xi_1,add_Xi_1,p_2,n_2,Xi_2,add_Xi_2);

new_w = (T*w')';
Wf = diag(1./new_w);

P = Wf*T*Wc;
