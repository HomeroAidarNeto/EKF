function [yout,xmodelo,x0default,H,u,MM,MMy] = modelocin(x0,ti,t,xinf0)
//Modelo cinético de Schultz 2020 com parâmetros da Nicole - simulados dia 03/03/2023

    exec('lib_ekf/ctegamma.sci');
    exec('lib_ekf/eqdif.sci');
    exec('lib_ekf/vel.sci');

    //Parâmetros do modelo - 
    kcat = 1.2723e+07; // 1/min
    kcatl = 6.6458e+09;
    Kmh = 1.6832e-05;
    Kmt = 3.9631e-06;
    kh = 773.65;
    kt = 2.2287e+05;
    Ki = 1.6721e-05;
    Km = 2.1837;
    Kmgal = 2.7413e-06
    ke=   4.7821e-06 // 1/min 
    ctevel= [kcat,kcatl,Kmh,Kmt,kh,kt,Ki,Km,Kmgal,ke]
    
    //Entradas do sistema (u)
    CAe  = 0.5     //mol L
    Q    = 0.000    // L/min testando batelada n alimentada
    u    = [CAe,Q]
    
    //Estado inicial (mol/L) (x0)
    if xinf0 == 0 then
        Lac = 0.31 //se não temos concentrações reais, usa o valor do artigo
    else
        Lac  = xinf0(1,1) 
    end
    Glu  = 0
    Gal  = 0
    Glb  = 0  //GOS 2
    Tri  = 0  //GOS 3
    Trig = 0  //GOS 3
    Tet  = 0  //GOS 4
    Tetg = 0  //GOS 4
    Ets = [ 2.00E-06  
            2.00E-07
            6.00e-8  //4.00E-07/4
            4.00E-07
            4.00E-07
            2.00E-06
            2.00e-8  //2.00E-07
            9.00E-08
            2.00E-07]
    if xinf0 == 0 then
        Et = 4.00E-07 //se não temos concentrações reais, usa um valor qualquer
    else
        Et  = Ets(xinf0(1,$)) // ultimo elemento de xinf0 tem numero do ensaio
    end
    V0   = 1 //L
 
    x0default   = [Lac,Glu,Gal,Glb,Tri,Trig,Tet,Tetg,Et,V0]'
    if x0 == 0 then
        x0 = x0default
    end
    xmodelo = ode(x0,ti,t,list(eqdif,u,ctevel))'
    //['di', 'gli', 'gal', 'gos3','gos4']
    H = [1 0 0 1 0 0 0 0 0 0;
         0 1 0 0 0 0 0 0 0 0;
         0 0 1 0 0 0 0 0 0 0;
         0 0 0 0 1 1 0 0 0 0;
         0 0 0 0 0 0 1 1 0 0;
         1 0 0 0 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0 0];
         
    yout = xmodelo*H'
    
    MM = [2*180-18, 180, 180, 2*180-18,3*180-2*18,3*180-2*18,4*180-3*18,4*180-3*18,9e+4,1] //errado?
    MMy= [MM(4),MM(2),MM(3),(MM(5)+MM(6))/2,(MM(7)+MM(8))/2,MM(1),MM(4)]
endfunction
