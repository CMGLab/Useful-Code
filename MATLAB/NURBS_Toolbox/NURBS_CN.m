function [] = NURBS_CN(d,P)

if d ==2
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
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    hold on
    axis equal
    
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
return