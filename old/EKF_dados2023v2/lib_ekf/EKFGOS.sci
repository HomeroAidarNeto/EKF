function  EKFGOS(C,cname,texp,t,xinf0)
    
    exec('lib_ekf/jacobiano.sci');
    exec('lib_ekf/ctegamma.sci');
    exec('lib_ekf/eqdif.sci');
    exec('lib_ekf/vel.sci');
    
    selecaoy= [4,6]         
    selecao = [1,3,5,6] // seleciona quais espécies mostrar nos gráficos
    especies= ['Lac','Glu','Gal','Glb','Tri','Trig','Tet','Tetg','Et']

    [yteo,xmod,x0,H,u,MM,MMy] = modelocin(0,0,t,xinf0) 
    
     yexp = C //saidas medidas, Xinferido PLS
    
    // plot conc molar
    scf(1)
    plot(t,xmod(:,1:$-1))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(especies)
    title('Molar x tempo - todas as espécies')
    
    scf(2)
    plot(t,yteo(:,selecaoy));plot(texp,yexp(:,selecaoy),'o');plot(texp,xinf0(:,selecaoy),'x')
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(cname(selecaoy))
    title('Molar x tempo - seleção')
    
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
    title('Concentrações teo e exp x tempo') 
       
    // plot conc massica
    xmodM= xmod.*repmat(MM,length(t),1) //converte de mol/L pra g/L
    yexpM = yexp.*repmat(MMy,length(texp),1)
//    scf(4)
//    plot(t,xmodM(:,1:$-1))
//    xlabel("t (min)", "fontsize", 4);
//    ylabel("C (g/L)", "fontsize", 4);
//    legend(['Lac','Glu','Gal','Glb','Tri','Trig','Tet','Tetg','Et'])
//    title('Mássica x tempo - todas as espécies')
//    
//    scf(5)
//    plot(t,xmodM(:,selecao))
//    xlabel("t (min)", "fontsize", 4);
//    ylabel("C (g/L)", "fontsize", 4);
//    legend(especies(selecao))
//    title('Mássica x tempo - seleção')
  
//    scf(4)
//    plot(texp,yexpM,'o')
//    scf(5)
//    plot(texp,yexpM(:,selecao),'o')
    
    
// Variáveis auxiliares
    mr = size(xmod,'r') // número de pontos de dados
    mexp = size(yexp,'r') // número de pontos experimentais
    nx = size(x0,'r') //número de variáveis de estado  
    nu = length(u) // número de variáveis de entrada
    ny = size(yexp,'c') // número de variáveis de saída

//Variáveis para o EKF    
    sigy = [0.0150127   0.0063762   0.0019683   0.0044827   0.0003461   0.0150189   0.0000625] //Obtido rodando main-multi com mesmo pretrat e dados - RMSECV pro numero otimo de regressores (mol/L)
    sigx    = ones(1,nx)*0.01
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
    if t(i)>= texp(ix) then// tenho uma amostra
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
        if t(i)>= texp(ix) then// tenho uma amostra
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
        plot(t(1:i),Ypost(1:i,selecaoy),'-')
        plot(t(1:i),Ypostm(1:i,selecaoy),'--')
        plot(t(1:i),YpostM(1:i,selecaoy),'--')
    
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
                
        
        
//        scf(4)

//        plot(t(1:i),Xpost(1:i,1:ny).*repmat(MM(1:$-1),length(t(1:i)),1),'-r')
//        plot(t(1:i),Xpostm(1:i,1:ny).*repmat(MM(1:$-1),length(t(1:i)),1),'--')
//        plot(t(1:i),XpostM(1:i,1:ny).*repmat(MM(1:$-1),length(t(1:i)),1),'--')
//        
//        scf(5)
//        plot(t(1:i),Xpost(1:i,selecao).*repmat(MM(selecao),length(t(1:i)),1),'-r')
//        plot(t(1:i),Xpostm(1:i,selecao).*repmat(MM(selecao),length(t(1:i)),1),'--')
//        plot(t(1:i),XpostM(1:i,selecao).*repmat(MM(selecao),length(t(1:i)),1),'--')  
    end // end aqui plota tudo em tempo real, mais lento
endfunction 
