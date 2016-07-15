%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Transform_Single_to_Multi
%
% Input:  d = number of spatial dimensions
%         P = array of NURBS control points (single-indexed)
%         w = array of NURBS weights (single-indexed)
%         n_1 = number of basis functions in direction 1
%         n_2 = number of basis functions in direction 2
%
% Output: P_multi = array of NURBS control points
%                   (multi-indexed)
%         w_multi = array of NURBS weights
%                   (multi-indexed)
%
% Purpose: Map control points to multi-index
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [P_multi,w_multi] = Transform_Single_to_Multi(d,P,w,n_1,n_2)

P_multi = zeros(n_1,n_2,d);
for i = 1:d
    P_multi(:,:,i) = reshape(P(:,i),[n_2,n_1])';
end

w_multi = reshape(w,[n_2,n_1])';