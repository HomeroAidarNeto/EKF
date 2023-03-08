clc; clear ; for i = 1:20; close ; end;

exec('lib_multi/func_multi_ILS_new.sci');
exec('lib_ekf/EKFGOS.sci');
exec('lib_ekf/modelocin.sci');
    exec('lib_multi/diffmeu.sci');
    exec('lib_multi/alisar.sci');
    exec('lib_multi/pls.sci');
    exec('lib_multi/spa_clean.sci');
    exec('lib_multi/spa_model.sci');
    exec('lib_multi/pls_model.sci');
    exec('lib_multi/pcr_model.sci');
    exec('lib_multi/zscore.sci');
    exec('lib_multi/loess.sci');
    exec('lib_multi/polyfit.sci');
    exec('lib_multi/func_pretreatment.sci');
    exec('lib_multi/func_analysis.sci');
    exec('lib_multi/sgolay_filt.sci');

kmax = 10; // está relacionado com os componentes que estão contidos na amostra (não colocar muito )

cname = ['di', 'gli', 'gal', 'gos3','gos4'];

x0=[]
absor0 = []
for i=[3:3] // 1 e 8 estão ruins
    if i<10 then
        arqent = '0'+string(i)+'.dat'
    else
        arqent = string(i)+'.dat'
    end
    xlendo = fscanfMat('x'+arqent);
    absorlendo = fscanfMat('abs'+arqent); 
    x0=[x0;xlendo]
    absor0=[absor0;absorlendo(2:$,:)]
end
absor0 = [absorlendo(1,:);absor0]

klb = x0\absor0(2:$,250:470)
t   =[0:1:size(x0,'r')-1]'
yout = modelocin(t)
Absout  = yout*klb //tem elementos maiores que 1 e negativos??
[Yp,Ytp,Par_norm] = pls_model(absor0(2:$,250:470),x0,kmax,Absout)
yexp = Yp
scf(1); plot(yexp,yout);

[Beta,T]=pls(Absout,x0,kmax);
yexp2 = Absout*Beta;
pause
clf(1)
EKFGOS(yexp,cname)


