function xdot = eqdif(t,x,u,ctevel) // Balanço de massa
    exec('lib_ekf/vel.sci');
    CAe  = u(1)      // entradas do sistema 
    Q    = u(2)
    
    Lac  = x(1)
    Glu  = x(2)
    Gal  = x(3)
    Glb  = x(4)
    Tri  = x(5)
    Trig = x(6)
    Tet  = x(7)
    Tetg = x(8)
    Et   = x(9)
    v    = x(10)
    
    r  = vel(x,ctevel) 
    rLac = r(1)
    rGlu = r(2)
    rGal = r(3)
    rGlb = r(4)
    rTri = r(5)
    rTrig= r(6)
    rTet = r(7)
    rTetg= r(8)
    rEt  = r(9)
    
    if t>80 then //Corta alimentação em t=80
        Q=0
    end
    
    Lacdot  = (CAe*Q + rLac* v - Lac* Q)./v 
    Gludot  = (  0   + rGlu* v - Glu* Q)./v
    Galdot  = (  0   + rGal* v - Gal* Q)./v
    Glbdot  = (  0   + rGlb* v - Glb* Q)./v
    Tridot  = (  0   + rTri* v - Tri* Q)./v
    Trigdot = (  0   + rTrig*v - Trig*Q)./v
    Tetdot  = (  0   + rTet* v - Tet* Q)./v
    Tetgdot = (  0   + rTetg*v - Tetg*Q)./v
    Etdot   = (  0   + rEt*  v - Et*  Q)./v  
    Vdot    = Q
  
    xdot = [Lacdot Gludot Galdot Glbdot Tridot Trigdot Tetdot Tetgdot Etdot Vdot]
endfunction
