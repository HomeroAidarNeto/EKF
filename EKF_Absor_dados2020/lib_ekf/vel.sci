//Velocidades de reação
function r = vel(concs,ctevel)
    exec('lib_ekf/ctegamma.sci');
    Lac  = concs(1) 
    Glu  = concs(2)
    Gal  = concs(3)
    Glb  = concs(4)
    Tri  = concs(5)
    Trig = concs(6)
    Tet  = concs(7)
    Tetg = concs(8)
    Et   = concs(9)
    [Gamma, a, b, E]= ctegamma(concs,ctevel)
   
    kcat = ctevel(1)
    kcatl= ctevel(2)
    Kmh  = ctevel(3)
    Kmt  = ctevel(4)
    kh   = ctevel(5)
    kt   = ctevel(6)
    Ki   = ctevel(7)
    Km   = ctevel(8)
    Kmgal= ctevel(9)
    ke   = ctevel(10)
    
    rLac = E*(kh/Kmh *Tri - Lac*(kcat/Km + Gamma*kt/Kmt)) 
    rGlu = E*(kcat*Lac/Km)
    rGal = E*(Gamma*(kcatl - kt/Kmgal*Gal)+kh/Kmh*Glb)
    rGlb = E*(Gamma*kt*(Gal/Kmgal-Glb/Kmt)+kh/Kmh*(Trig-Glb))
    rTri = E*(Gamma*kt/Kmt*(Lac-Tri)+kh/Kmh*(Tet-Tri))
    rTrig= E*(Gamma*kt/Kmt*(Glb-Trig)+kh/Kmh*(Tetg-Trig))
    rTet = E*(Gamma*kt/Kmt*Tri-kh/Kmh*Tet)
    rTetg= E*(Gamma*kt/Kmt*Trig-kh/Kmh*Tetg)
    rEt  = -ke*Et
    
    r = [rLac,rGlu,rGal,rGlb,rTri,rTrig,rTet,rTetg,rEt]
endfunction
