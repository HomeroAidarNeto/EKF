function [absor,lambda,x,ifig] = pretrat(absor0,x0,pretreat,analysis)
    
    lambda0=absor0(1,:); // first row is the wavelength (or wavenumber)
    absor0=absor0(2:$,:);
    absor = absor0;
    lambda=lambda0;
    x = x0;
    ifig = 1;

    [absor,lambda,x,ifig]=func_pretreatment(pretreat,absor,lambda,x,ifig)

    [absor,lambda,x,ifig]=func_analysis(analysis,absor,lambda,x,ifig)
endfunction
