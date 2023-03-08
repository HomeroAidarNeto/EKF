//Constantes a,b e gamma
function [Gamma, a, b, E] = ctegamma(concs,ctevel)
    Lac  = concs(1)
    Glu  = concs(2)
    Gal  = concs(3)
    Glb  = concs(4)
    Tri  = concs(5)
    Trig = concs(6)
    Tet  = concs(7)
    Tetg = concs(8)
    Et   = concs(9)
    
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
    
a    = Glb+Tri+Trig+Tet+Tetg
b    = Lac+Glb+Tri+Trig
Gamma= (kcat/Km * Lac + kh/Kmh * a)/(kcatl+kt*(Gal/Kmgal + b/Kmt))
E       = Et/(1+Gal/Ki+Lac/Km+a/Km+Gamma*(1+Gal/Kmgal+b/Kmt))

endfunction
