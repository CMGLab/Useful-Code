function [N] = NURBS_1D_T(xi, i, p,n, Xi,w)

N = zeros(n,1);

for j = 1:n
    N(j) = deBoor(xi,i(j),p,Xi);   
end


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