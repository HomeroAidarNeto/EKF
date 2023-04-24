
//   ----------------------------------------------------------- //
//    MULTI_ILS.M            
//
//   ultima atualizaçao 2020
//   Marcelo P. A. Ribeiro
//   ----------------------------------------------------------- //

// Este programa AJUSTA e VALIDA - UTILIZANDO VALIDAÇAO CRUZADA - multicalibraçao
// utilizando o metodo de Regressao Linear Multipla (multiple linear regression,
// MLR ou, tambem chamada de inverse Least Square, ILS).
// Dois métodos para a reduçao da singularidade da matriz A'A podem ser
// utilizadas.
//
// A entrada ocorre em dois arquivos: 1o arquivo contém a matriz com
// concentrações, onde a linha é a amostra e a coluna é o componente; o
// 2o arquivo é uma matriz onde a primeira coluna é o número do comprimento,
// depois as linhas indicam amostras e colunas o valor de absorbância em
// cada comprimento de onda.
//
// Neste programa apenas um conjunto de dados deve ser inserido.
// Deste conjunto serao tirados os dados de ajuste e validaçao.
//

////==========================================================================================================////
////                                   INICIO DE INSERÇAO DE DADOS                                            ////
////==========================================================================================================////

clc; clear ; for i = 1:20; close ; end;

exec('lib_multi/func_multi_ILS_new_conjuntos.sci');
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

//// Método
Selecao = 1; // 1 = PLS; 2 = SPA; 3=PCR ;4 = Lolimot
// Para o SPA é necessário dar o primeiro comprimento de onda a ser
// utilizado, lini.
optkini = 2; // 0=> lini = lambda(1); 1=> lini = dado abaixo; 2=> otimiza lini.
lini = 0; // [nm]. Só tem sentido se optkini = 1

//// número máximo de regressores para avaliar na  Validacao Cruzada
kmax = 15; // está relacionado com os componentes que estão contidos na amostra (não colocar muito )

//------------------------------------------------------------------
//Número de Analitos (colunas no arquivo de concentrações)
// -----------------------------------------------------------------
nc = 7;
//------------------------------------------------------------------
//Nome dos Analitos e unidade de concentração
// -----------------------------------------------------------------
cname = ['di', 'gli', 'gal', 'gos3','gos4','lac','glb'];
unid = 'M';

// --------------------------------------------------------------
//Carregando os dados
// --------------------------------------------------------------
// x = matriz de concentracoes [nc x nd]
// absor = matriz com as absorbancias, sendo a primeira linha o comprimento
// de onda [nd+1 x nl]

//Matriz 3D com todos os dados
aux = zeros(1,211*7*9);
xtot = matrix(aux,211,7,9);
aux = zeros(1,(211+1)*911*9);
absortot = matrix(aux,(211+1),911,9);
for i=[1:9] 
    if i<10 then
        arqent = '0'+string(i)+'.dat'
    else
        arqent = string(i)+'.dat'
    end
    xtot(:,:,i) = fscanfMat('x'+arqent);
    absortot(:,:,i) = fscanfMat('abs'+arqent); 
end

// Fração de dados aleatórios para testes
frac_test =.1 ; //em geral de 0.0 a 0.4

// Optimization of the model complexity
// 1) Validation data constant.
// Fraction of remaining data used for validation (after test data removal).
frac_val =.20 ; //usually 0.05-0.4
OptimModel = {'Val',frac_val}
// 2) CrossValidation ('k-fold'): kpart = number of partitions. 
//    if kpart = -1 then kpart = nd (leave one out - LOOCV)
kpart = 4;
OptimModel = {'kfold',kpart}
rho = 2 // Pro akaike
//OptimModel = {'IC','AIC',rho}
// 3) Information criteria 
// Akaike: {'IC','AIC',rho); Baysian: {'IC','BIC'}; Khinchin: {'IC','LILC',rho}; 
// Final prediction: {'IC','FPE'}; Strutctural risk min: {'IC','SRM',rho1,rho2}
// rho's are tuned parameters.

// Set of pretreatment 
// Moving Average: {'MA',radius,Losing points = 1 or 2, plot=1 or 0}
// LOESS: {'Loess',alpha = [0.2-0.9], order = [1,2,...],plot=1 or 0}
// Savitzky and Golay Filter: {'SG',radius,Losing points = 1 or 2,poly order = integer, der order = integer, plot=1 or 0}   
// Derivative: {'Deriv',order,plot=1 or 0}
// Cutting regions: {'Cut',lower bound,upper bound,plot=1 or 0}
// Cutting maxAbs: {'CutAbs',maxAbs, action: warning = 0 or cutting = 1 ,plot=1 or 0}
// Control Baseline based on noninformation region: {'BLCtr',ini_lamb,final_lamb,Abs_value, plot=1 or 0}
// example: 
//pretreat = list({'SG',6,1,3,1,1},{'Cut',270,470,1});// controle de baseline todos
pretreat = list({'SG',6,1,3,1,1},{'Cut',300,400,1});// controle de baseline todos

// Set of Analysis
// Lambert-Beer: {'LB'}
// Principal Component Analysis (PCA): {'PCA'} 
analysis = list({'LB'},{'PCA'});
//analysis = list({'PCA'})
//analysis = []

// outlier analysis
outlier = 1; // 0 = no; 1 = yes.

// Gravar Xval para os diferentes número de regressores
gravarXval = 0; //1 = grava, 0 = não grava (=1 pode dar erro)

for index = [1:9]
    intervalo=[]
    for i = [1:9]
        if i ~= index
            intervalo = [intervalo,i];
        end
    end
   dadosteste = {xtot(:,:,index),absortot(:,:,index)}
   
   x0=[]
   absor0=[]
   for i = intervalo
       x0=[x0;xtot(:,:,i)]
       absor0=[absor0;absortot(2:$,:,i)]
   end
   absor0 = [absortot(1,:,1);absor0]
   
for i = 1:20; close ; end;
    [Xp,Xpval,Par_norm,absor,x] = func_multi_ILS_new_conjuntos(Selecao,optkini,lini,kmax,nc,cname,unid,x0,absor0,frac_test,dadosteste,OptimModel,pretreat,analysis,outlier,gravarXval,index)
end



