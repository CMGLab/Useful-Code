function [N] = NURBS_2D(xi_1, xi_2, i_1,i_2,p_1,p_2,n_1,n_2,Xi_1,Xi_2,w)
N_1 = zeros(n_1,1);
N_2 = zeros(n_2,1);

for j = 1:n_1
    N_1(j) = deBoor(xi_1,j,p_1,Xi_1);
end

for j = 1:n_2
    N_2(j) = deBoor(xi_2,j,p_2,Xi_2);
end

[N_2g,N_1g] = meshgrid(N_2,N_1);

if size(N_2g)~=size(w)
    N_1g = N_1g';
    N_2g = N_2g';
end

R = N_1(i_1)*N_2(i_2)*w(i_1,i_2)/(sum(sum(N_1g.*N_2g.*w)));

N  = R;


return

function [N] = deBoor(xi, i, p, Xi)
if p==0
    
    if xi==Xi(i) && xi==Xi(i+1)
        N = 0;
    elseif xi>=Xi(i) && xi<=Xi(i+1)
        N = 1;
    else
        N = 0;
    end
    
else
    if Xi(i) ==Xi(i+p)
        f = 0;
    else
        f = (xi - Xi(i))/(Xi(i+p)-Xi(i));
    end
    
    if Xi(i+1) == Xi(i+p+1)
        g = 0;
    else
        g = (Xi(i+p+1)-xi)/(Xi(i+p+1)-Xi(i+1));
    end
    
    N = f*deBoor(xi,i,p-1,Xi) + g*deBoor(xi,i+1,p-1,Xi);
    
end

return