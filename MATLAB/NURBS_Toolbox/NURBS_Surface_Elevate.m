function [new_p_1,new_p_2,new_n_1,new_n_2,new_Xi_1,new_Xi_2,new_P,new_w] =...
NURBS_Surface_Elevate(d,dir,p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w)
% Do all the calculation in projective space: That is multiply the control
% points by their weights:
for i = 1:d
    P(:,:,i) = P(:,:,i).*w;
end
P(:,:,d+1) = w;

if dir ==1
    nel = length(unique(Xi_1))-1;    
    ctr1 = 1;
    ctr2 = 1;
    
    for e = 1:nel
        [Q(ctr1:ctr1+p_1+1,:,:),new_p_1] = elevateDegree(P(ctr2:ctr2+p_1,:,:),p_1,dir);
        ctr1 = ctr1+p_1+1;
        ctr2 = ctr2+p_1;
    end
        
    XiLoc = unique(Xi_1);
    nXi = size(Q,1)+new_p_1+1;
    % Go in and change the knot vector.
    % Initialize the desired final knot vector:
    new_Xi_1 = zeros(1,nXi);
    
    % The first knot will have multiplicity of 4, so don't change the first
    % four entries. The last knot will also have multiplicity of 4, so
    % assign the last four entries to 1
    new_Xi_1(nXi-new_p_1:nXi) = ones(1,new_p_1+1);
    
    % The middle knots will all have multiplicity of p_1, so loop through and
    % assign the rest of the knots
    ctr = new_p_1+2;
    for iKV = 2:length(XiLoc)-1
        for j = 0:new_p_1-1
            new_Xi_1(ctr+j)   = XiLoc(iKV);
        end
        ctr = ctr+new_p_1;
    end

    new_n_1 = n_1+nel;
    
    new_p_2 = p_2;
        new_n_2 = n_2;
            new_Xi_2 = Xi_2;


elseif dir ==2
    nel = length(unique(Xi_2))-1;    
    ctr1 = 1;
    ctr2 = 1;
    
    for e = 1:nel
        [Q(:,ctr1:ctr1+p_2+1,:),new_p_2] = elevateDegree(P(:,ctr2:ctr2+p_2,:),p_2,dir);
        ctr1 = ctr1+p_2+1;
        ctr2 = ctr2+p_2;
    end
        
    XiLoc = unique(Xi_2);
    nXi = size(Q,2)+new_p_2+1;
    % Go in and change the knot vector.
    % Initialize the desired final knot vector:
    new_Xi_2 = zeros(1,nXi);
    
    % The first knot will have multiplicity of 4, so don't change the first
    % four entries. The last knot will also have multiplicity of 4, so
    % assign the last four entries to 1
    new_Xi_2(nXi-new_p_2:nXi) = ones(1,new_p_2+1);
    
    % The middle knots will all have multiplicity of p_1, so loop through and
    % assign the rest of the knots
    ctr = new_p_2+2;
    for iKV = 2:length(XiLoc)-1
        for j = 0:new_p_2-1
            new_Xi_2(ctr+j)   = XiLoc(iKV);
        end
        ctr = ctr+new_p_2;
    end
    new_n_2 = n_2+nel;

        new_p_1 = p_1;
        new_n_1 = n_1;
            new_Xi_1 = Xi_1;

end

new_w = Q(:,:,d+1);
for i = 1:d
    new_P(:,:,i) = Q(:,:,i)./new_w;
end

return
    
function [Pe,pe] = elevateDegree(P,p,dir)

% This function performs degree elevation on a NURBS curve of degree p.
% It takes as input the control points of the existing nurbs curve, P and its
% knot vector and outputs the control points and knot vector if a degree p+1
% NURBS curve.

if dir == 1
    n = size(P,dir);
    
    ne = n+1;
    pe = p+1;
    
    Pe(1,:,:) = P(1,:,:);
    Pe(ne,:,:) = P(n,:,:);
    
    for i = 2:n
        Pe(i,:,:) = (i-1)/(p+1)*P(i-1,:,:) + (1 - (i-1)/(p+1))*P(i,:,:);
    end
    
elseif dir ==2
    n = size(P,dir);
    
    ne = n+1;
    pe = p+1;
    
    Pe(:,1,:) = P(:,1,:);
    Pe(:,ne,:) = P(:,n,:);
    
    for i = 2:n
        Pe(:,i,:) = (i-1)/(p+1)*P(:,i-1,:) + (1 - (i-1)/(p+1))*P(:,i,:);
    end
end
return