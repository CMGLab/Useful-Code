function [] = NURBS_Curve(d,p,n,Xi,P,w,varargin)

if nargin == 6
    CN = true;
    span = [0 max(Xi)];
    LW = 1;

elseif nargin == 7
    CN = varargin{1};
    span = [0 max(Xi)];
    LW = 1;
    
elseif nargin == 8
    CN = varargin{1};
    span = varargin{2};
    LW = 1;
elseif nargin == 9
    CN = varargin{1};
    span = varargin{2};
    LW = varargin{3};
end

h = .002;
xi  = span(1):h:span(2);

C = zeros(length(xi),d);

for j = 1:length(xi)
    for i = 1:n
        C(j,:) = C(j,:) + NURBS_1D(xi(j), i, p,n,Xi,w)*P(i,:);
    end    
end

if d ==2
    
    plot(C(:,1),C(:,2),'b','LineWidth',LW)
    hold on
    if CN
      scatter(P(:,1),P(:,2),'or','filled')
      plot(P(:,1),P(:,2),'k')
    end
    xlabel('X')
    ylabel('Y')

elseif d == 3
    figure
    plot3(C(:,1),C(:,2),C(:,3),'b','LineWidth',LW)
    hold on
    if CN
    scatter3(P(:,1),P(:,2),P(:,3),'or','filled')
    plot3(P(:,1),P(:,2),P(:,3),'k')
    end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
else
    disp('Error: This function only supports curves of dimension 2 or 3.')
end
return