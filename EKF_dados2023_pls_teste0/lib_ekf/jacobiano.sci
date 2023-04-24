function J = jacobiano(x0,ti,tf,xpri,nx,xinf0)
    J = zeros(nx,nx)
    hx = 1e-7
    for j=1:nx // calculando a matriz A = dx/dx0
        x0j=x0
        x0j(j)=x0(j)+hx // incremento no n=j-esimo estado
        [lixo,xhxj] = modelocin(x0j,t(i-1),t(i),xinf0)
        xhxj = xhxj'
        J(:,j) = (xhxj-xpri)/hx
    end    
endfunction
