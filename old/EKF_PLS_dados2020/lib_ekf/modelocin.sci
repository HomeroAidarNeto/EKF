function [yout,xmodelo,x0default,H,u,MM,MMy] = modelocin(x0,ti,t,xinf0)
//Modelo cinético com parametros de Schultz 2020 (Guilhermina)

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
//  Et   = 2.42e-7/2 artigo
    Et   = 4e-8 
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
         0 0 0 0 0 0 1 1 0 0];
         
    yout = xmodelo*H'
    
    MM = [2*180-18, 180, 180, 2*180-18,3*180-2*18,3*180-2*18,4*180-3*18,4*180-3*18,9e+4,1]
    MMy= [MM(4),MM(2),MM(3),(MM(5)+MM(6))/2,(MM(7)+MM(8))/2]
endfunction
