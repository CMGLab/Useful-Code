function [IEN] = IEN_1D(p,n,Xi)

% IEN_1D
% Written by Christopher Coley
% Last modified 24 Sep 14
%
% Computes the IEN array for a one-dimensional B-spline basis
%
% Output:
%   IEN - array mapping local basis function/element numbers to global
%   basis function numbers
%
% Input:
%   p - polynomial degree of the basis function to be evaluated; scalar
%   n - total number of basis functions; scalar
%   Xi - knot vector; row vector

l = p+1;
e = 1;

while l < n+1
    for a = 1:p+1
        IEN(a,e) = (l+a) - (p+1);
    end
    
    l = l+1;
    
    while Xi(l+1) == Xi(l) && l < n+1
        l = l+1;
    end
    
    if l < n+1
        e = e+1;
    end
end