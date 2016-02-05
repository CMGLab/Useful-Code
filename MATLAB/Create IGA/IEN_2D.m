function [IEN] = IEN_2D(p_1,p_2,n_1,n_2,Xi_1,Xi_2)

% IEN_2D
% Written by Christopher Coley
% Last modified 24 Sep 14
%
% Computes the IEN array for a tw0-dimensional B-spline basis
%
% Output:
%   IEN - array mapping local basis function/element numbers to global
%   basis function numbers
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

IEN_1 = IEN_1D(p_1,n_1,Xi_1);
IEN_2 = IEN_1D(p_2,n_2,Xi_2);

n_el_1 = size(IEN_1,2);
n_el_2 = size(IEN_2,2);

IEN = zeros(p_1*p_2+p_1+p_2+1,n_el_1*n_el_2);

for e1 = 1:n_el_1
    for a1 = 1:p_1+1
        
        i1 = IEN_1(a1,e1);
        
        for e2 = 1:n_el_2
            for a2 = 1:p_2+1
                
                i2 = IEN_2(a2,e2);
                e = (e1-1)*n_el_2 + e2;
                a = (a1-1)*(p_2+1) + a2;
                i = (i1-1)*n_2 + i2;
                
                IEN(a,e) = i;
                
            end
        end
    end
end