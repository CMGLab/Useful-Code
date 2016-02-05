function [ n_el,C_operators,IEN,P_b,w_b,P_e,w_e,n_el_1,n_el_2 ] = ExtractAndLocalize( p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w )
% ExtractAndLocalize
% Written by Christopher Coley and AJ Gemer
% Last modified 7 Oct 14
%
% Conducts Bezier extraction and element localization
%
% Output:
%   n_el - number of Bezier elements
%   C_operators - array storing element extraction operators
%   IEN - array mapping local basis function/element numbers to global
%   basis function numbers
%   P_b - array of localized Bezier points for each element
%   w_b - array of localized Bezier weights for each element
%   P_e - array of NURBS points for each element
%   w_b - array of NURBS weights for each element
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
%   P - Array of Bezier control points; 2 x (n_1 x n_2) array [xi;yi]
%   w - array storing NURBS weights; n_1 x n_2 matrix

P_e = [];
w_e = [];
n_loc = (p_1+1)*(p_2+1);
d = 2;

w = reshape(w',1,n_1*n_2);

[n_el,C_operators,IEN,n_el_1,n_el_2] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2);
    
[P_b,w_b] = Extract_Geometry(d,p_1,p_2,n_el,C_operators,IEN,P,w);
    
    for e = 1:n_el
        for a = 1:n_loc
            i = IEN(a,e);
            P_e = cat(2,P_e,P(:,i));
            w_e = cat(2,w_e,w(i)); 
        end
    end
end

