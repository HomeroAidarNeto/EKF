function [yout,xmodelo,x0default,H,u] = modelocin(x0,ti,t)

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
    Lac  = 0.315//dissacarideo
    Glu  = 5.268535e-003
    Gal  = 2.508022e-003
    Glb  = 0  //GOS 2
    Tri  = 2.616172e-003  //GOS 3
    Trig = 0  //GOS 3
    Tet  = 0  //GOS 4
    Tetg = 0  //GOS 4
//    Et   = 2.42e-7/2 
    Et   = 4e-8 
    V0   = 1 //L
    t0=0
    x0default   = [Lac,Glu,Gal,Glb,Tri,Trig,Tet,Tetg,Et,V0]'
    
    if  x0==0 then
        x0 = x0default
    end
    if ti == 0 then
        ti = t0
    end

    xmodelo = ode(x0,ti,t,list(eqdif,u,ctevel))'
    //['di', 'gli', 'gal', 'gos3','gos4']
    H = [1 0 0 1 0 0 0 0 0 0;
         0 1 0 0 0 0 0 0 0 0;
         0 0 1 0 0 0 0 0 0 0;
         0 0 0 0 1 1 0 0 0 0;
         0 0 0 0 0 0 1 1 0 0];
         
    yout = xmodelo*H'
endfunction
