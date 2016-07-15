%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Transform_Multi_to_Single
%
% Input:  d = number of spatial dimensions
%         P = array of NURBS control points (multi-indexed)
%         w = array of NURBS weights (multi-indexed)
%         n_1 = number of basis functions in direction 1
%         n_2 = number of basis functions in direction 2
%
% Output: P_single = array of NURBS control points
%                    (single-indexed)
%         w_single = array of NURBS weights
%                    (single-indexed)
%
% Purpose: Map control points from multi-index
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_single,w_single] = Transform_Multi_to_Single(d,P,w,n_1,n_2)

P_single = zeros(n_1*n_2,d);
for i = 1:d
    P_single(:,i) = reshape(transpose(P(:,:,i)),n_1*n_2,1);    
end

w_single = reshape(transpose(w),n_1*n_2,1);