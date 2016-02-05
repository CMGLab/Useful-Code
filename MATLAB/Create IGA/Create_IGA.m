% Create_IGA
% Written by Christopher Coley
% Last modified 10 Nov 15
%
% Generates a .iga file for simple geometries.
%
% Output:
%   N/A
%
% Input:
%   N/A

function [n_el,n_1,n_2] = Create_IGA(dirz,filename,p,nRefine)

full_filename = [dirz,filename,'.iga'];

%Define the geometry
if p == 2
    
    % Geometry for p = 2
    dim = 2;
    p = 2;
    n_1 = 3;
    n_2 = 3;
    Xi_1 = [0,0,0,1,1,1];
    Xi_2 = [0,0,0,1,1,1];
    P = [0,0.5,1,0,0.5,1,0,0.5,1;
        0,0,0,0.5,0.5,0.5,1,1,1];
    w = [1,1,1;
        1,1,1;
        1,1,1];
else
    
    % Geometry for p = 3
    dim = 2;
    p = 3;
    n_1 = 4;
    n_2 = 4;
    Xi_1 = [0,0,0,0,1,1,1,1];
    Xi_2 = [0,0,0,0,1,1,1,1]; 
    P = [0,1/3,2/3,1,0,1/3,2/3,1,0,1/3,2/3,1,0,1/3,2/3,1;
        0,0,0,0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,2/3,1,1,1,1];
    w = [1,1,1,1;...
        1,1,1,1;...
        1,1,1,1;...
        1,1,1,1];
    
end

%Refine the geometry
if nRefine > 0
    close all
    for i=1:nRefine
        add_Xi_1=[];
        add_Xi_2=[];
        for j=1:length(Xi_1)-1
            if Xi_1(j) ~= Xi_1(j+1)
                temp_xi1 = (Xi_1(j+1)-Xi_1(j))/2 + Xi_1(j);
                add_Xi_1 = [add_Xi_1 temp_xi1];
            end
        end
        for j=1:length(Xi_2)-1
            if Xi_2(j) ~= Xi_2(j+1)
                temp_xi2 = (Xi_2(j+1)-Xi_2(j))/2 + Xi_2(j);
                add_Xi_2 = [add_Xi_2 temp_xi2];
            end
        end
        
        %Reshape variables to match expected input for NURBS_Surface_Refine
        P = ShapePoints(P,n_1,n_2,dim);
        
        %Refine surface
        [n_1,n_2,Xi_1,Xi_2,P,w] = NURBS_Surface_Refine(dim,add_Xi_1,add_Xi_2,p,p,n_1,n_2,Xi_1,Xi_2,P,w);
        
        %Reshape variables to match expected input for Extract_Basis
        P = UnshapePoints(P,n_1,n_2,dim);
       
    end
end

%Perform Bezier extraction
[ n_el,C_operators,IEN,~,~,~,~,~,~ ] = ExtractAndLocalize( p,p,n_1,n_2,Xi_1,Xi_2,P,w );
IEN = IEN - 1;

%Create file and write data
fileID = fopen(full_filename,'wt');

if dim == 2
    fprintf(fileID,'%s\n','type plane');
else
    fprintf(fileID,'%s\n','type surface');
end

fprintf(fileID,'nodeN %d\n',numel(w));
fprintf(fileID,'elemN %d\n',n_el);

if dim == 2
    P = cat(1,P,zeros(1,numel(w)));
end
w = reshape(w',1,numel(w));
P = cat(1,P,w);
fprintf(fileID,'node %.15f %.15f %.15f %.15f\n',P);

n_loc = (p+1)^2;

for i = 1:n_el
    fprintf(fileID,'belem %d %d %d\n',length(IEN(:,1)),p,p);
    fprintf(fileID,'%d ',IEN(:,i));
    fprintf(fileID,'\n');
    for j = 1:n_loc
        fprintf(fileID,'%.15f ',C_operators(j,:,i));
        fprintf(fileID,'\n');
    end
end

fclose(fileID);