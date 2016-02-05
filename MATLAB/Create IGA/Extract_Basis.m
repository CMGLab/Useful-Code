function [n_el,C_operators,IEN,n_el_1,n_el_2] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2)

% Extract_Basis
% Written by Christopher Coley
% Last modified 14 Oct 14
%
% Constructs the element extraction operators and corresponding IEN array
% for a two-dimensional B-spline basis
%
% Output:
%   n_el - number of Bezier elements
%   C_operators - array storing element extraction operators
%   IEN - array mapping local basis function/element numbers to global
%   basis function numbers
%   n_el_1 - number of Bezier elements in first direction
%   n_el_2 - number of Bezier elements in second direction
%
% Input:
%   p_1 - polynomial degree of the basis function to be evaluated in the
%   first dimension; scalar
%   p_2 - polynomial degree of the basis function to be evaluated in the
%   second dimension; scalar
%   n_1 - total number of basis functions in the first dimension; scalar
%   n_2 - total number of basis functions in the second dimension; scalar
%   Xi_1 - knot vector in the first dimension; row vector
%   Xi_2 - knot vector in the second dimension; row vector

[n_el_1,C_operators_1,IEN_1] = Extract_Basis_1D(p_1,n_1,Xi_1);
[n_el_2,C_operators_2,IEN_2] = Extract_Basis_1D(p_2,n_2,Xi_2);

n_el = n_el_1*n_el_2;

for e2 = 1:n_el_2
    for e1 = 1:n_el_1
        e = (e2-1)*n_el_1 + e1;
        C_operators(:,:,e) = kron(C_operators_2(:,1+(e2-1)*(p_2+1):p_2+1+(e2-1)*(p_2+1)),C_operators_1(:,1+(e1-1)*(p_1+1):p_1+1+(e1-1)*(p_1+1)));
    end
end

IEN = IEN_2D(p_1,p_2,n_1,n_2,Xi_1,Xi_2);