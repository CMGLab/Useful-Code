function [Pnew] = UnshapePoints(P,n_1,n_2,d)

% UnhapePoints
% Written by Christopher Coley and AJ Gemer
% Last modified 15 Oct 14
%
% Reshapes points from our Luke's NURBS_Surface_Refine code's format 
% to be compatible with our code
%
% Output:
%   Pnew - Array of reshaped points
%
% Input:
%   P - array of points
%   p_1 - polynomial order of the element in the first dimension; scalar
%   p_2 - polynomial order of the element in the second dimension; scalar
%   d - spatial dimensions; scalar

Pnew = zeros(d,n_1*n_2);
for i = 1:d
    Pnew(i,:) = reshape(P(:,:,i)',[n_1*n_2,1])';
end