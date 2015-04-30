function [P,new_n_1,new_Xi_1,new_n_2,new_Xi_2] = Prolongation_Operator_2D(p_1,n_1,Xi_1,add_Xi_1,p_2,n_2,Xi_2,add_Xi_2)

% Prolongation_Operator_2D
% Written by Christopher Coley
% Last modified 1 Apr 15
%
% Creates the prolongation operator for 2D using B-splines
%
% Output:
%   P - 2D prolongation operator
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
% Input:
%   p_1 - polynomial degree of the original and refined B-spline curves in
%   first dimension;  scalar
%   n_1 - number of basis functions for the original B-spline space in first
%   dimension; scalar
%   Xi_1 - univariate knot vector for the original B-spline space in first
%   dimension; row vector
%   add_Xi_1 - knots to be added in first dimension; row vector
%   p_2 - polynomial degree of the original and refined B-spline curves in
%   second dimension;  scalar
%   n_2 - number of basis functions for the original B-spline space in second
%   dimension; scalar
%   Xi_2 - univariate knot vector for the original B-spline space in second
%   dimension; row vector
%   add_Xi_2 - knots to be added in second dimension; row vector

[P1,new_n_1,new_Xi_1] = Prolongation_Operator_1D(p_1,n_1,Xi_1,add_Xi_1);
[P2,new_n_2,new_Xi_2] = Prolongation_Operator_1D(p_2,n_2,Xi_2,add_Xi_2);

P = kron(P1,P2);