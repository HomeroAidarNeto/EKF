//   ----------------------------------------------------------- //
//    INFER_ILS.M
//
//   ultima atualizaçao 2020
//   Marcelo P. A. Ribeiro
//   ----------------------------------------------------------- //

// Este programa INFERE a concentração a partir de Regressao Linear Multipla
// (multiple linear regression, MLR). E utiliza os dados inferidos para
// a entrada do EKF, para inferir o estado do reator.
// Três métodos para a reduçao da singularidade da matriz A'A podem ser
// utilizadas. Neste código, utiliza-se o PLS.
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
Selecao = 1; // 1 = PLS; 2 = SPA; 3=PCR; 4 = Lolimot
optkini = 2; // 0=> lini = lambda(1); 1=> lini = dado abaixo; 2=> otimiza lini.
lini = 0; // [nm]. Só tem sentido se optkini = 1

//// número de regressores para inferência
kinf = [4 1 3 1 1 4 1]; //kinf está relacionado com os componentes que estão contidos na amostra (não colocar muito )
//------------------------------------------------------------------
//Número de Analitos (colunas no arquivo de concentrações)
// -----------------------------------------------------------------
nc = 7;
//------------------------------------------------------------------
//Nome dos Analitos e unidade de concentração
// -----------------------------------------------------------------
cname = ['di', 'gli', 'gal', 'gos3','gos4','lac', 'glb'];
unid = 'M';

// --------------------------------------------------------------
// Carregando os dados de Ajuste
// --------------------------------------------------------------
// x = matriz de concentracoes [nc x nd]
// absor = matriz com as absorbancias, sendo a primeira linha o comprimento
// de onda [nd+1 x nl]
x0=[]
absor0 = []
for i=[7:9] 
    if i<10 then
        arqent = '0'+string(i)+'.dat'
    else
        arqent = string(i)+'.dat'
    end
    xlendo = fscanfMat('x'+arqent);
    absorlendo = fscanfMat('abs'+arqent); 
    x0=[x0;xlendo]
    absor0=[absor0;absorlendo(2:$,:)]
    nexp = size(xlendo,'r')
end
absor0 = [absorlendo(1,:);absor0]
// --------------------------------------------------------------
// Carregando os dados para Inferência
// --------------------------------------------------------------
// xinf = matriz de concentracoes [nc x nd] (OPCIONAL)
// absorinf = matriz com as absorbancias, sendo a primeira linha o comprimento
// de onda [nd+1 x nl]
xinf0=[]
absorinf0 = []
for i=[7] 
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

// Escolha do pre tratamento 
// Savitzky and Golay Filter: {'SG',radius,Losing points = 1 or 2,poly order = inte
//pretreat = list({'SG',6,1,3,1,1},{'Cut',270,470,1});// controle de baseline todos
pretreat = list({'SG',6,1,3,1,1},{'Cut',270,350,1});// controle de baseline todos
pretreatinf = pretreat;

// Escolha da análise
analysis = [];
analysisinf = [];

// Outlier Analysis
leverage = 0; // 1 = faz análise; 0 = não faz análise

//Definindo variáveis para a função de inferência
texp0   =[0:1:nexp]' 
t      = [0:1:nexp]' //tempo do modelo começa em 0 
texp  = texp0([1:1:nexp]) //tempo experimental pega texp0 de 5 em 5 minutos, começa em 1 pq?

absorinf0 = absorinf0([1 [1:1:nexp]+1],:) //1° linha = lambdas, depois pega absorbancias em cada tempo experimental, como se fosse uma sonda medindo de n em n pontos, pro xinf0 é a memsa coisa
xinf0 = xinf0([1:1:nexp],:) //isto supõe que as concentrações e absorbancias corrspondem ao tempo, 600 medidas demoraram 600 minutos, como arrumamos isso?

[Xp,Xinf,Par_norm,absor,x,Xinf_conc] = func_infer_ILS(Selecao,optkini,lini,kinf,nc,cname,unid,x0,absor0,xinf0,absorinf0,pretreat,pretreatinf,analysis,analysisinf,leverage);

pause
for index = 1:100; close; end;
xinf0(1,$+1) = i // xinf0 carrega o numero do ensaio, util p definir cond inicial de enzima real
EKFGOS(Xinf_conc,cname,texp,t,xinf0)  
