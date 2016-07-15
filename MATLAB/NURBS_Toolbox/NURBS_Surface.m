function [] = NURBS_Surface(d,p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,h,varargin)

if nargin == 9
    CN = true;
else
    CN = varargin{1};
end

xi_1  = 0:h:max(Xi_1);
xi_2  = 0:h:max(Xi_2);

S = zeros(length(xi_1),length(xi_2),d);



for j_1 = 1:length(xi_1)
    for j_2  = 1:length(xi_2)
        for i_1 = 1:n_1
            for i_2 = 1:n_2
                S(j_1,j_2,:) =  S(j_1,j_2,:) + ...
                    NURBS_2D(xi_1(j_1), xi_2(j_2), i_1,i_2,...
                    p_1,p_2,n_1,n_2,Xi_1,Xi_2,w)*...
                    P(i_1,i_2,:);
            end
        end
    end
end
if d ==2
    figure
    surf(S(:,:,1),S(:,:,2),zeros(size(S(:,:,1))),'EdgeAlpha',0,'FaceColor','g')
    xlabel('X')
    ylabel('Y')
    view(0,90)
    hold on
    axis equal
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            scatter(P(i,j,1),P(i,j,2),'or','filled')
        end
    end
    for i = 1:size(P,1)-1
        for j = 1:size(P,2)-1
            plot([P(i,j,1),P(i+1,j,1)],[P(i,j,2),P(i+1,j,2)],'k')
            plot([P(i,j,1),P(i,j+1,1)],[P(i,j,2),P(i,j+1,2)],'k')
            
        end
        j = size(P,2);
        plot([P(i,j,1),P(i+1,j,1)],[P(i,j,2),P(i+1,j,2)],'k')
        
    end
elseif d == 3
    %figure
    surf(S(:,:,1),S(:,:,2),S(:,:,3),'EdgeAlpha',0,'FaceColor','g')%,'FaceLighting','phong')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    hold on
    axis equal
    if CN
        for i = 1:size(P,1)
            for j = 1:size(P,2)
                scatter3(P(i,j,1),P(i,j,2),P(i,j,3),'or','filled')
            end
        end
        
        for i = 1:size(P,1)-1
            for j = 1:size(P,2)-1
                plot3([P(i,j,1),P(i+1,j,1)],[P(i,j,2),P(i+1,j,2)],...
                    [P(i,j,3),P(i+1,j,3)],'k')
                plot3([P(i,j,1),P(i,j+1,1)],[P(i,j,2),P(i,j+1,2)],...
                    [P(i,j,3),P(i,j+1,3)],'k')
                
            end
            j = size(P,2);
            plot3([P(i,j,1),P(i+1,j,1)],[P(i,j,2),P(i+1,j,2)],...
                [P(i,j,3),P(i+1,j,3)],'k')
            
        end
        
        i  = size(P,1);
        
        
        for j = 1:size(P,2)-1
            plot3([P(i,j,1),P(i,j+1,1)],[P(i,j,2),P(i,j+1,2)],...
                [P(i,j,3),P(i,j+1,3)],'k')
        end
        
    end
else
    disp('Error: This function only supports curves of dimension 2 or 3.')
end
return