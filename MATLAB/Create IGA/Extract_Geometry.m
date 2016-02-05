function [P_b,w_b] = Extract_Geometry(d,p_1,p_2,n_el,C_operators,IEN,P,w)

% Extract_Geometry
% Written by Christopher Coley
% Last modified 24 Sep 14
%
% Computes the Bezier control points and weights corresponding to a NURBS
% surface
%
% Output:
%   P_b - array storing the Bezier control points for the NURBS surface
%   w_b - array storing the Bezier weights for the NURBS surface
%
% Input:
%   d - spatial dimension of the physical space; scalar
%   p_1 - polynomial degree of the basis function to be evaluated in the
%   first dimension; scalar
%   p_2 - polynomial degree of the basis function to be evaluated in the
%   second dimension; scalar
%   n_el - number of Bezier elements; scalar
%   C_operators - array storing element extraction operators; n_el x n_el
%   x n_el array (?)
%   IEN - array mapping local basis function/element numbers to global
%   basis function numbers; n_el x n_el array (?)
%   P - array of control points; d x (n_1 x n_2) array [xi;yi;zi]
%   w - array storing NURBS weights; n_1 x n_2 array

%multiply by weights and concatenate weights to control points
w = reshape(w',1,[]);
Pw = zeros(d+1,length(w));
for j = 1:d
    Pw(j,:) = P(j,:).*w;
end
Pw(d+1,:) = w;

for e = 1:n_el
    for j = 1:d+1
        for a = 1:(p_1+1)*(p_2+1)
            i = IEN(a,e);
            P_e(a,1) = Pw(j,i);
        end
        C_e = C_operators(:,:,e);
        P_bw(j,:,e) = (C_e'*P_e)';
    end
end

%divide out the weights and extract new control points and weights
for e = 1:n_el
    w_b(1,:,e) = P_bw(d+1,:,e);
    for j = 1:d
        P_b(j,:,e) = P_bw(j,:,e)./w_b(:,:,e);
    end
end

