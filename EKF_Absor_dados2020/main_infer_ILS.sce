//   ----------------------------------------------------------- //
//    INFER_ILS.M
//
//   ultima atualizaçao 2020
//   Marcelo P. A. Ribeiro
//   ----------------------------------------------------------- //

// Este programa INFERE a concentração a partir de Regressao Linear Multipla
// (multiple linear regression, MLR).
// Três métodos para a reduçao da singularidade da matriz A'A podem ser
// utilizadas.
//
// Para o ajuste:
// A entrada ocorre em dois arquivos: 1o arquivo contém a matriz com
// concentrações, onde a linha é a amostra e a coluna é o componente; o
// 2o arquivo é uma matriz onde a primeira linha é o valor do comprimento de onda,
// as linhas abaixo indicam amostras e colunas o valor de absorbância em
// cada comprimento de onda.
//
// Para a inferência:
// Um arquivo de varreduras, onde a primeira linha indica o comprimento de
// onda. Um arquivo de concentração poderá ser dado também, se quiser
// avaliar o modelo em dados de teste externo.
//
////==========================================================================================================////
////                                   INICIO DE INSERÇAO DE DADOS                                                        ////
////==========================================================================================================////

clc; clear ; for i = 1:100; close ; end;

exec('lib_multi/func_infer_ILS.sci');
exec('lib_multi/func_multi_ILS_new.sci');
exec('lib_ekf/EKFGOS.sci');
exec('lib_ekf/EKFGOS2.sci');
exec('lib_ekf/EKFGOS3.sci');
exec('lib_ekf/modelocin.sci');
exec('lib_ekf/pretrat.sci');
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


//// Método
Selecao = 1; // 1 = PLS; 2 = SPA; 3=PCR
// Para o SPA é necessário dar o primeiro comprimento de onda a ser
// utilizado, lini.
optkini = 2; // 0=> lini = lambda(1); 1=> lini = dado abaixo; 2=> otimiza lini.
lini = 0; // [nm]. Só tem sentido se optkini = 1

//// número de regressores para inferência
//kinf = [2 2 2 3 2]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) UV
kinf = [2 3 3 3 3];//[3 3 3 1 3]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) UV der
//kinf = [4 4 4 3 3]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) NIRA der
//kinf = [2 2 2 3 1]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) NIRA der
//kinf = [10 10 10 10 10]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) UV 2019

//kinf = [4 4 4 4 4]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) UV der todos
//kinf = [5 5 5 5 4]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) UV der todos 2
//kinf = [4 4 4 4 4]; // está relacionado com os componentes que estão contidos na amostra (não colocar muito ) UV  todos 2
kinf = [2 2 2 2 6]; // UV  der 2019 
kinf = [4 4 5 4 4]; // UV  todos 2 SG
kinf = [5 5 5 5 5]; // UV  der todos 2 SG
//------------------------------------------------------------------
//Número de Analitos (colunas no arquivo de concentrações)
// -----------------------------------------------------------------
nc = 5;
//------------------------------------------------------------------
//Nome dos Analitos e unidade de concentração
// -----------------------------------------------------------------
cname = ['di', 'gli', 'gal', 'gos3','gos4'];
unid = 'M';

// --------------------------------------------------------------
// Carregando os dados de Ajuste
// --------------------------------------------------------------
// x = matriz de concentracoes [nc x nd]
// absor = matriz com as absorbancias, sendo a primeira linha o comprimento
// de onda [nd+1 x nl]


//x0 = fscanfMat('X_Multi_UV2019.txt');     
//absor0 = fscanfMat('Abs_Multi_UV2019.txt'); 
// x0 = load('C_Dados_puros.txt');     // Estes arquivos foram pre-tratados pelo
// absor0 = load('Dados_puros.txt'); // Dados para teste
//
//x0 = fscanfMat('X_Multi_UV2020.txt');     
//absor0 = fscanfMat('Abs_Multi_UV2020.txt'); 

//x0 = fscanfMat('X_Multi_NIRA2020.txt');     
//absor0 = fscanfMat('Abs_Multi_NIRA2020.txt'); 
//
//x0(1,:)=[]
//absor0(2,:)=[] 

//Script original
//x0 = fscanfMat('X_Multi_UV2020.txt');     
//absor0 = fscanfMat('Abs_Multi_UV2020.txt'); 
//x00 = fscanfMat('X_Multi_UV2019.txt');     
//absor00 = fscanfMat('Abs_Multi_UV2019.txt');
//
//x0 = [x0 ; x00]
//absor0 = [absor0 ; absor00(2:$,:)] 

x0=[]
absor0 = []
for i=[2] // 1 e 8 estão ruins
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

// --------------------------------------------------------------
// Carregando os dados para Inferência
// --------------------------------------------------------------
// xinf = matriz de concentracoes [nc x nd] (OPCIONAL)
// absorinf = matriz com as absorbancias, sendo a primeira linha o comprimento
// de onda [nd+1 x nl]
//xinf0 =[];
////xinf0 = fscanfMat('X_Multi_UV2020.txt');     
////absorinf0 = fscanfMat('Abs_Multi_UV2020.txt'); 
//
//xinf0 = x0;     
//absorinf0 = absor0; 

xinf0=[]
absorinf0 = []
for i=[6] // 1 e 8 estão ruins
    if i<10 then
        arqent = '0'+string(i)+'.dat'
    else
        arqent = string(i)+'.dat'
    end
    xlendo = fscanfMat('x'+arqent);
    absorlendo = fscanfMat('abs'+arqent); 
    xinf0=[xinf0;xlendo]
    absorinf0=[absorinf0;absorlendo(2:$,:)]
end
absorinf0 = [absorlendo(1,:);absorinf0]

//xinf0 = fscanfMat('X_Multi_NIRA2020.txt');     
//absorinf0 = fscanfMat('Abs_Multi_NIRA2020.txt'); 

// Set of pretreatment 
// Moving Average: {'MA',radius,Losing points = 1 or 2, plot=1 or 0}
// LOESS: {'Loess',alpha = [0.2-0.9], order = [1,2,...],plot=1 or 0}
// Savitzky and Golay Filter: {'SG',radius,Losing points = 1 or 2,poly order = integer, der order = integer, plot=1 or 0}   
// Derivative: {'Deriv',order,plot=1 or 0}
// Cutting regions: {'Cut',lower bound,upper bound,plot=1 or 0}
// Cutting maxAbs: {'CutAbs',maxAbs, action: warning = 0 or cuting = 1 ,plot=1 or 0}
// example: 
    //pretreat = list({'MA',3,1,1}, {'Deriv',1,1},{'MA',3,1,1},{'Loess',0.2,1,1},{'Cut',240,350, 355,360,1});

pretreat = list({'Cut',270,550,1},{'MA',6,1,1});
pretreat = list({'Cut',270,550,1},{'MA',6,1,1}, {'Deriv',1,1},{'MA',6,1,1});
pretreat = list({'MA',6,1,1}, {'Deriv',1,1},{'MA',6,1,1},{'Cut',270,470,1}); // der todos 2
pretreat = list({'SG',6,1,3,1,1},{'Cut',270,470,1}); // der todos 2
pretreat = list({'SG',9,1,3,0,1},{'Cut',270,470,1});//  todos 2

pretreat = list({'SG',9,1,2,1,1},{'Cut',270,470,1});// der todos 2
pretreat = list({'SG',6,1,3,1,1},{'Cut',270,470,1});// controle de baseline todos
//pretreat = list({'MA',6,1,1}, {'Deriv',1,1},{'MA',6,2,1});
//pretreat = list({'MA',6,1,1},{'Cut',270,470,1}); // todos 
pretreatinf = pretreat;

// Set of Analysis
// Lambert-Beer: {'LB'}
// Principal Component Analysis (PCA): {'PCA'} 
// analysis = list({'LB'},{'PCA'});
// analysisinf = [];

analysis = list({'LB'},{'PCA'});
analysis = [];
analysisinf = []; //list({'LB'},{'PCA'});

// Outlier Analysis
// Leverage
leverage = 0; // 1 = faz análise; 0 = não faz análise



nd = size(x0,'r')
//xone = [ones(nd,1) x0];
xone =  x0;

indice = find(absor0(1,:)>=270 &  absor0(1,:)<=470     )
klb = xone\absor0(2:$,indice) //ajuste com erro Eabs
// se fizermos klb com termo indep, klb - nc+1 x nlamb, H nc+1 x nlamb+1, x nx+1 x 1
Abssim = xone*klb
Eabs = absor0(2:$,indice)- Abssim
gdl = size(Abssim,'c') * ( size(Abssim,'r') - size(klb,'r') )//npontos - nparametros
S2 = sum(Eabs.^2)/(gdl )
sigy = sqrt(S2)
texp0   =[1:1:600]'

t      = [0:1:600]'
texp  = texp0([1:5:600])
xinf0 = xinf0([1:5:600],:)
absorinf0 = absorinf0([1 [1:5:600]+1],:)

[yout] = modelocin(0,0,texp) //Adicionar ruido da medida aqui ou na absout?
[yout] = modelocin(0,0,t)
nd = size(yout,'r')

//Absout  = [ones(nd,1) yout]*klb
Absout  =  yout*klb
//Absout = Absout + rand(Absout,'n')*stdev(Eabs) //tem elementos maiores que 1 e negativos??

Absout = [absor0(1,indice); Absout]

//xinf0 = yout;    
//absorinf0 = Absout; 


//absorinf0 = [absor0(1,indice); absor0(texp,indice)] 
//xinf0 = x0(texp,:)

//absorinf0($+1,:) = absor0($,indice)//faz com que  ultimo do absorinf seja iigual ultimo do absor0
//xinf0($+1,:) = x0($,:)


//Abssim = [absor0(1,indice); Abssim]

//[Xp,Xinf,Par_norm,absor,x,Xinf_conc] = func_infer_ILS(Selecao,optkini,lini,kinf,nc,cname,unid,x0,absor0,xinf0,absorinf0,pretreat,pretreatinf,analysis,analysisinf,leverage);


//clf(1); scf(1); plot(Xinf_conc,yout);
//pause
for i = 1:100; close; end;
//EKFGOS(Xinf_conc,cname,texp,t) 

//Absor = absor0(2:$,indice)
//Absor = Absor(texp,:)
//EKFGOS2(Absor,cname,texp,t,klb) 

[Absor,lambda,x,ifig] = pretrat(absor0,x0,pretreat,analysis)
klb2 = x\Absor //ajuste com um conjunto de dados (2)
    scf(5); plot(yout(2:$,:),x,[0,.35],[0,.35])
    xlabel("C modelo (mol/L)", "fontsize", 4);
    ylabel("C pretrat HPLC (mol/L)", "fontsize", 4);
    title('Conc modelo vs Conc do conjunto de dados HPLC')

absor0 = [] //Coleta absorbancias "medidas" pra utilizar no EKF
for i=[4] // 1 e 8 estão ruins
    if i<10 then
        arqent = '0'+string(i)+'.dat'
    else
        arqent = string(i)+'.dat'
    end
    absorlendo = fscanfMat('abs'+arqent); 
    absor0=[absor0;absorlendo(2:$,:)]
end
absor0 = [absorlendo(1,:);absor0]

[Absorexp,lambda,x,ifig] = pretrat(absor0,x0,pretreat,analysis) //tem que pretratar essas absorbâncias "medidas"?

pause
for i = 1:100; close; end;
EKFGOS3(Absorexp,cname,texp,t,klb2,sigy,lambda,x0) 
