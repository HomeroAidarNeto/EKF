function  EKFGOS2(Absor,cname,texp,t,klb)
    
    exec('lib_ekf/jacobiano.sci');
    exec('lib_ekf/ctegamma.sci');
    exec('lib_ekf/eqdif.sci');
    exec('lib_ekf/vel.sci');
    
    //Parâmetros do modelo
    kcat  = 3.526e+7 // 1/min
    kcatl = 9.943e+9 // 1/min
    Kmh   = 3.345e-6 // mol/L
    Kmt   = 1.758e-6 // mol/L
    kh    = 2.876e+3 // 1/min
    kt    = 1.439e+4 // 1/min
    Ki    = 9.405e-6 // mol/L
    Km    = 2.586e-2 // mol/L
    Kmgal = 2.433e-6 // mol/L
    ke    = 1.651e-4 // 1/min
    ctevel= [kcat,kcatl,Kmh,Kmt,kh,kt,Ki,Km,Kmgal,ke]
    
    //Entradas do sistema (u)
    CAe  = 0.5     //mol L
    Q    = 0.000    // L/min testando batelada n alimentada
    u    = [CAe,Q]
    
    //Estado inicial (mol/L) (x0)
    Lac  = 0.41//dissacarideo
    Glu  = 0
    Gal  = 0
    Glb  = 0  //GOS 2
    Tri  = 0  //GOS 3
    Trig = 0  //GOS 3
    Tet  = 0  //GOS 4
    Tetg = 0  //GOS 4
//    Et   = 2.42e-7/2 
    Et   = 4e-8 
    V0   = 1 //L
    x0   = [Lac,Glu,Gal,Glb,Tri,Trig,Tet,Tetg,Et,V0]' 
        
    //Vetor de massa molar para 9 espécies, ultimo elemento é pra acertar a dimensão da matriz ao usar repmat.
    MM = [2*180-18, 180, 180, 2*180-18,3*180-2*18,3*180-2*18,4*180-3*18,4*180-3*18,9e+4,1]
    MMy= [MM(4),MM(2),MM(3),(MM(5)+MM(6))/2,(MM(7)+MM(8))/2]
    
    selecao = [2,3,4] // seleciona quais espécies mostrar nos gráficos
    especies= ['Lac','Glu','Gal','Glb','Tri','Trig','Tet','Tetg','Et']

    [yteo,xmod,x0,H,u] = modelocin(t) 
    xmodM= xmod.*repmat(MM,length(t),1) //converte de mol/L pra g/L
    // plot conc molar
    scf(1)
    plot(t,xmod(:,1:$-1))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(['Lac','Glu','Gal','Glb','Tri','Trig','Tet','Tetg','Et'])
    title('Molar x tempo - todas as espécies')
    
    scf(2)
    plot(t,xmod(:,selecao))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(especies(selecao))
    title('Molar x tempo - seleção')
    
    scf(3)
    plot(t,xmod(:,$))
    xlabel("t (min)", "fontsize", 4);
    ylabel("V (L)", "fontsize", 4);
    legend(['V'])
    title('Volume no reator x tempo')
    
    // plot conc massica
    scf(4)
    plot(t,xmodM(:,1:$-1))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (g/L)", "fontsize", 4);
    legend(['Lac','Glu','Gal','Glb','Tri','Trig','Tet','Tetg','Et'])
    title('Mássica x tempo - todas as espécies')
    
    scf(5)
    plot(t,xmodM(:,selecao))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (g/L)", "fontsize", 4);
    legend(especies(selecao))
    title('Mássica x tempo - seleção')
    
    scf(7)
    plot(t,xmod(:,$-1))
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (g/L)", "fontsize", 4);
    legend("Enzima")
    title('Molar x tempo - seleção')
 
    C = Absor*klb'
    yexp = C //saidas medidas
   
    yexpM = yexp.*repmat(MMy,length(texp),1)
    scf(4)
    plot(texp,yexpM,'o')
//    scf(5)
//    plot(texp,yexpM(:,selecao),'o')
    
    mr = size(xmod,'r') // número de pontos de dados
    mexp = size(yexp,'r') // número de pontos experimentais
    nx = size(x0,'r') //número de variáveis de estado  
    nu = length(u) // número de variáveis de entrada
    ny = size(yexp,'c') // número de variáveis de saída
    
    sigy = [0.0329816   0.0187083   0.0117486   0.0134983   0.0014765] //Obtido rodando main-multi com mesmo pretrat e dados - RMSECV pro numero otimo de regressores (mol/L)
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

    scf(6); plot(texp,yexp,'o',t,yteo)
    legend(cname)
    scf(6); plot(t,yteo)
    xlabel("t (min)", "fontsize", 4);
    ylabel("C (mol/L)", "fontsize", 4);
    legend(cname)
    title('Concentrações teo e exp x tempo')
    
    Xpost = zeros(xmod) // declarando variável
    XpostM = zeros(xmod)// declarando variável
    Xpostm = zeros(xmod)// declarando variável
    i=1 //tempo zero t(1)
    ix=1 // 1a amostra 
    
    //xpri = yexp(1,:)' // chutando o próprio dado experimental
    //xpri = [yexp(1,:)'; V0+rand(y(1,nx),'normal')*sig]
    xpri = x0 // chutando o valor real

    K = Ppri*H'*pinv(H*Ppri*H'+V*R*V')
    xpost = xpri + K*(yexp(ix,:)'-H*xpri) // utilizando como chute o 1o dado exp
    ind   = find(xpost<0) // impede estados negativos
    xpost(ind) = 0
    ix=2 // indice do proximo dado experimental
    Ppost = (eye(nx,nx)-K*H)*Ppri // Calc a covariancia do x posterior
    EPx = sqrt(diag(Ppost)) // Erro padrão do estado posterior
    // guardando dados
    Xpost(i,:) = xpost' 
    XpostM(i,:) = xpost'+ 2*EPx' // intervalo de ~95% de confiança
    Xpostm(i,:) = xpost'- 2*EPx'// intervalo de ~95% de confiança
    
    for i=2:mr //Loop Kalman
        // Time  update
        dt = t(i)-t(i-1)
        x0 = xpost
        xpri = ode(x0,t(i-1),t(i),list(eqdif,u,ctevel)) // x prior - estimativa do modelo
        A = jacobiano(x0,t(i-1),t(i),xpri,nx)
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
        plot(t(1:i),Xpost(1:i,selecao),'-')
        plot(t(1:i),Xpostm(1:i,selecao),'--')
        plot(t(1:i),XpostM(1:i,selecao),'--')
    
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
