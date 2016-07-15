function [new_n,new_Xi,new_P,new_w] = NURBS_Curve_Refine(d,add_Xi,p,n,Xi,P,w)

% Determine what new_Xi should look like:
new_Xi = sort([Xi, add_Xi]);

% Do all the calculation in projective space: That is multiply the control
% points by their weights:
for i = 1:d
    P(:,i) = P(:,i).*w;
end
P(:,d+1) = w;

% Knot Interstion:
% We know what we want our final knot vector to look like, so we'll do knot
% insertion until our KV matches our KVF.
for i = 1:length(new_Xi)
    % Check to see if the current knot in KV matches the current knot in
    % KVF. If it does, great, move on to the next knot, if not, calculate
    % new control points and insert KVF(i) into KV(i).
    if new_Xi(i) == Xi(i)
        continue
    else
        
        % Correct for the algorithm indexing by zero
        ki = i-1;
        
        % Calculate the new control points
        % Loop through indexes from k-p+1 to k and calculate the new control
        % points
        for j = ki-p+1:ki
            a(j) = (new_Xi(i)-Xi(j))/(Xi(j+p)-Xi(j));
            Q(j,:) = (1-a(j))*P(j-1,:) + a(j)*P(j,:);
        end
        
        % Move control points below ki-1 down by 1 space to
        % make room for the p new points;
        P(ki+1:length(P)+1,:) = P(ki:length(P),:);
        % Insert the new control points
        P(ki-p+1:ki,:) = Q(ki-p+1:ki,:);
        
        % move the knots to the right of the knot to be inserted over one
        % index
        Xi(i+1:end+1) = Xi(i:end);
        
        % Insert the knot to the knot vector.
        Xi(i) = new_Xi(i);
    end
end
new_w = P(:,d+1);
for i = 1:d
    new_P(:,i) = P(:,i)./new_w;
end

new_n = size(new_P,1);
return 