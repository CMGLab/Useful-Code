function [P,new_n,new_Xi] = Prolongation_Operator_1D(p,n,Xi,add_Xi)

% Prolongation_Operator
% Written by Christopher Coley
% Last modified 2 Apr 15
%
% Creates the prolongation operator for 1D using B-splines
%
% Output:
%   P - 1D prolongation operator
%   new_n - number of basis functions for the refined B-spline space
%   new_Xi - univariate knot vector for the refined B-spline space
%
% Input:
%   p - polynomial degree of the original and refined B-spline curves; scalar
%   n - number of basis functions for the original B-spline space; scalar
%   Xi - univariate knot vector for the original B-spline space; row vector
%   add_Xi - knots to be added; row vector

% initialize variables
P = eye(n,n);
new_n = n;

% insert one knot at a time and find the prolongation operator for that
% knot
% then multiply the current prolongation operator with the system
% prolongation operator 

for i = 1:length(add_Xi)
    
    % initialize the current prolongation operator
    new_n = new_n+1;
    P1 = zeros(new_n,n);
    
    % find k such that xi exists in [Xi(k), Xi(k+1))
    [TF,k]=ismember(0,add_Xi(i)<Xi,'legacy');
    
    for j = 1:n+i
        if j <= k-p
            P1(j,j) = 1;
        elseif j > k-p && j <= k

            % calculate alpha and insert into P
            a = (add_Xi(i)-Xi(j))/(Xi(j+p)-Xi(j));
            
            P1(j,j-1) = 1-a;
            P1(j,j) = a;
        else
            P1(j,j-1) = 1;
        end
    end
    
    new_Xi = [Xi(1:k) add_Xi(i) Xi(k+1:end)];
    
    Xi = new_Xi;
    
    P = P1*P;
end