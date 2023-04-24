function  EKFGOS(C,cname,texp,t,xinf0)
    
    exec('lib_ekf/jacobiano.sci');
    exec('lib_ekf/ctegamma.sci');
    exec('lib_ekf/eqdif.sci');
    exec('lib_ekf/vel.sci');
             
    selecao = [6,2,4,5] // seleciona quais espécies mostrar (y, cname)
    especies= ['Lac','Glu','Gal','Glb','Tri','Trig','Tet','Tetg','Et'] //(x, estado)

    [yteo,xmod,x0,H,u,MM,MMy] = modelocin(0,0,t,xinf0) 
    cenzima0 = xmod(1,9)
    yexp = C //saidas medidas, Xinferido PLS
    
    // plot conc molar
    scf(1)
    plot(t,xmod(:,1:$-1))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(especies)
    title('Molar x tempo - todas as espécies')
     
    scf(2); plot(texp,yexp(:,selecao),'o')
    legend(cname(selecao))
    scf(2); plot(t,yteo(:,selecao))
    legend(cname(selecao))
    scf(2); plot(texp,xinf0(:,selecao),'x')
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(cname(selecao))
    title(string('Concentrações medidas pelo softsensor (o), de teste (x), do modelo cinético (-), do EKF (- com linha pontilhada para o intervalo de confiança) vs tempo. Concentração inicial de enzima nos dados de teste = ')+string(cenzima0)+string(' mol/L') )
    
    scf(3)
    plot(t,xmod(:,$))
    xlabel("t (min)", "fontsize", 4);
    ylabel("V (L)", "fontsize", 4);
    legend(['V'])
    title('Volume no reator x tempo')
    
    scf(7)
    plot(t,xmod(:,$-1))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend("Enzima")
    title('Molar x tempo - seleção')
    
//    scf(6); plot(texp,yexp,'o',t,yteo,'x',t,xinf0)
    scf(6); plot(texp,yexp,'o',t,yteo,texp,xinf0(:,1:$-1),'x')
    legend(cname)
    scf(6); plot(t,yteo)
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(cname)
    title(string('Concentrações medidas pelo softsensor (o), de teste (x), do modelo cinético (-), do EKF (- com linha pontilhada para o intervalo de confiança) vs tempo. Concentração inicial de enzima nos dados de teste = ')+string(cenzima0)+string(' mol/L') )
           
// Variáveis auxiliares
    mr = size(xmod,'r') // número de pontos de dados
    mexp = size(yexp,'r') // número de pontos experimentais
    nx = size(x0,'r') //número de variáveis de estado  
    nu = length(u) // número de variáveis de entrada
    ny = size(yexp,'c') // número de variáveis de saída

//Variáveis para o EKF    
    sigy = [0.0083066   0.0018963   0.0012777   0.0027216   0.0000868   0.0035774   0.0000469] //Obtido rodando main-multi com mesmo pretrat e dados - RMSECV pro numero otimo de regressores (mol/L)
    sigx    = ones(1,nx)*0.001
    sigx(1,9) = 1e-8
    
    Ppri = diag([sigx.^2])
    sigw = 0.001// desvio padrão do ruído no processo
    
    q = eye(nx,nx)*sigw.^2 // covariancia do ruido do processo
    q(9,9) = (sigw*1e-7)^2
    R = diag(sigy.^2) // covarianacia do ruido das medidas
    
    A = zeros(nx,nx) //dx/dx0 - só declarando variável
    W = eye(nx,nx) // dx/dw - dependencia linear direta
    V = eye(ny,ny) // dy/dv - dependencia linear direta
      
    Xpost = zeros(xmod) // declarando variável
    XpostM = zeros(xmod)// declarando variável
    Xpostm = zeros(xmod)// declarando variável
    i=1 //tempo zero t(1)
    ix=1 // 1a amostra 
    
    //xpri = yexp(1,:)' // chutando o próprio dado experimental
    //xpri = [yexp(1,:)'; V0+rand(y(1,nx),'normal')*sig]
    xpri = x0 // chutando o valor real

    K = Ppri*H'*pinv(H*Ppri*H'+V*R*V')
    if ix <= size(texp,'r') & t(i)>= texp(ix) then// tenho uma amostra
        K = Ppri*H'*pinv(H*Ppri*H'+V*R*V')
        xpost = xpri + K*(yexp(ix,:)'-H*xpri) // x posterior
        Ppost = (eye(nx,nx)-K*H)*Ppri // covariancia do x posterior
        ix=ix+1 // indice próxima amostra
    else  // não tenho amostra
        xpost = xpri 
        Ppost = Ppri
    end
    ind   = find(xpost<0) // impede estados negativos
    xpost(ind) = 0
    EPx = sqrt(diag(Ppost)) // Erro padrão do estado posterior
    // guardando dados
    Xpost(i,:) = xpost' 
    XpostM(i,:) = xpost'+ 2*EPx' // intervalo de ~95% de confiança
    Xpostm(i,:) = xpost'- 2*EPx'// intervalo de ~95% de confiança
    
    Ypost(i,:) = Xpost(i,:)*H'
    YpostM(i,:) = XpostM(i,:)*H'
    Ypostm(i,:) = Xpostm(i,:)*H'
    
    for i=2:mr //Loop Kalman
        // Time  update
        dt = t(i)-t(i-1)
        x0 = xpost
        [lixo,xpri] = modelocin(x0,t(i-1),t(i),xinf0)
        xpri = xpri'
        A = jacobiano(x0,t(i-1),t(i),xpri,nx,xinf0)
        Ppri = A*Ppost*A'+W*q*W'// covariancia do x prior
        //Measurement Update
        if ix <= size(texp,'r')&t(i)>= texp(ix) then// tenho uma amostra
            K = Ppri*H'*pinv(H*Ppri*H'+V*R*V')
            xpost = xpri + K*(yexp(ix,:)'-H*xpri) // x posterior
            Ppost = (eye(nx,nx)-K*H)*Ppri // covariancia do x posterior
            ix=ix+1 // indice próxima amostra
        else  // não tenho amostra
            xpost = xpri 
            Ppost = Ppri
        end
        ind   = find(xpost<0)
        xpost(ind) = 0
       
        EPx = sqrt(diag(Ppost)) // Erro padrão
        // Guardando os dados 
        Xpost(i,:) = xpost'
        XpostM(i,:) = xpost'+ 2*EPx'
        Xpostm(i,:) = xpost'- 2*EPx'
        
        Ypost(i,:) = Xpost(i,:)*H'
        YpostM(i,:) = XpostM(i,:)*H'
        Ypostm(i,:) = Xpostm(i,:)*H'
//    end // end aqui plota tudo no final, mais rápido
        // plotando valores preditos pelo filtro
        scf(1)
        plot(t(1:i),Xpost(1:i,1:$-1),'-')
        plot(t(1:i),Xpostm(1:i,1:$-1),'--')
        plot(t(1:i),XpostM(1:i,1:$-1),'--')
        
        scf(2)
        plot(t(1:i),Ypost(1:i,selecao),'-')
        plot(t(1:i),Ypostm(1:i,selecao),'--')
        plot(t(1:i),YpostM(1:i,selecao),'--')
    
        scf(3)
        plot(t(1:i),Xpost(1:i,$),'-')
        plot(t(1:i),Xpostm(1:i,$),'--')
        plot(t(1:i),XpostM(1:i,$),'--')
        
        scf(7)
        plot(t(1:i),Xpost(1:i,$-1),'-')
        plot(t(1:i),Xpostm(1:i,$-1),'--')
        plot(t(1:i),XpostM(1:i,$-1),'--')
        
        scf(6)
        plot(t(1:i),Ypost(1:i,1:$),'-')
        plot(t(1:i),Ypostm(1:i,1:$),'--')
        plot(t(1:i),YpostM(1:i,1:$),'--')
                
    end // end aqui plota tudo em tempo real, mais lento
endfunction 
