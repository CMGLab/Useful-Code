function [n_el,C_operators,IEN] = Extract_Basis_1D(p,n,Xi)

% Extract_Basis_1D
% Written by Christopher Coley
% Last modified 24 Sep 14
%
% Constructs the element extraction operators and corresponding IEN array
% for a one-dimensional B-spline basis
%
% Output:
%   n_el - number of Bezier elements
%   C_operators - array storing element extraction operators
%   IEN - array mapping local basis function/element numbers to global
%   basis function numbers
%
% Input:
%   p - polynomial degree of the basis function to be evaluated; scalar
%   n - total number of basis functions; scalar
%   Xi - knot vector; row vector

a = p+1;
b = a+1;
n_el = 1;
C_current = eye(p+1);
C_operators= [];

while b<n+p+1
    C_next = eye(p+1);
    i = b;
    
    while b < n+p+1 && Xi(b+1) == Xi(b)
        b = b+1;
    end
    
    mult = b-i+1;
    
    if mult < p
        numer = Xi(b) - Xi(a);
        
        for j = p:-1:mult+1
            alphas(j-mult) = numer/(Xi(a+j)-Xi(a));
        end
        
        r = p-mult;
        
        for j = 1:r
            save = r-j+1;
            s = mult+j;
            
            for k = p+1:-1:s+1
                al = alphas(k-s);
                C_current(:,k) = al*C_current(:,k) + (1-al)*C_current(:,k-1);
            end
            
            if b < n+p+1
                C_next(save:j+save,save) = C_current(p-j+1:p+1,p+1);
            end
        end
    end
    
    C_operators = cat(2,C_operators,C_current);
    C_current = C_next;
    if b < n+p+1
        a = b;
        b = a+1;
        n_el = n_el+1;
    end
end

IEN = IEN_1D(p,n,Xi);