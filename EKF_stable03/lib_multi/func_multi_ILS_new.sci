function [Xp,Xpval,Par_norm,absor,x] = func_multi_ILS_new(Selecao,optkini,lini,kmax,nc,cname,unid,x0,absor0,frac_test,dadosteste,OptimModel,pretreat,analysis,outlier,gravarXval)

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
    flag = -1;
    ////==========================================================================================================////
    ////                                   INICIO DE INSERÇAO DE DADOS                                            ////
    ////==========================================================================================================////
    
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
    Par_norm = []
    //
    ////// Método
    //Selecao = 1; // 1 = PLS; 2 = SPA; 3=PCR
    //// Para o SPA é necessário dar o primeiro comprimento de onda a ser
    //// utilizado, lini.
    //optkini = 2; // 0=> lini = lambda(1); 1=> lini = dado abaixo; 2=> otimiza lini.
    //lini = 0; // [nm]. Só tem sentido se optkini = 1
    //
    ////// número máximo de regressores para avaliar na  Validacao Cruzada
    //kmax = 4; // está relacionado com os componentes que estão contidos na amostra (não colocar muito )
    //
    ////------------------------------------------------------------------
    ////Número de Analitos (colunas no arquivo de concentrações)
    //// -----------------------------------------------------------------
    //nc = 4;
    ////------------------------------------------------------------------
    ////Nome dos Analitos e unidade de concentração
    //// -----------------------------------------------------------------
    //cname = ['gos', 'glic', 'gal', 'lac'];
    //
    //
    //unid = 'g/L';
    //
    //// --------------------------------------------------------------
    ////Carregando os dados
    //// --------------------------------------------------------------
    //// x = matriz de concentracoes [nc x nd]
    //// absor = matriz com as absorbancias, sendo a primeira linha o comprimento
    //// de onda [nd+1 x nl]
    //
    //x0 = fscanfMat('XMultiCalExp.txt');     
    //absor0 = fscanfMat('AbsMultiCalExp_250-320.txt'); 
    //
    ////x0 = fscanfMat('Conc_aj-8-UVNIR.txt');     
    ////absor0 = fscanfMat('UV_aj_8.txt'); 
    ////absor0 = fscanfMat('NIR_aj_8.txt'); 
    //
    //// Fração de dados aleatórios para testes
    //frac_test =.0 ; //em geral de 0.0 a 0.4
    //
    //
    //// Set of pretreatment 
    //// Moving Average: {'MA',radius,Losing points = 1 or 2, plot=1 or 0}
    //// LOESS: {'Loess',alpha = [0.2-0.9], order = [1,2,...],plot=1 or 0}
    //// Derivative: {'Deriv',order,plot=1 or 0}
    //// Cutting regions: {'Cut',lower bound,upper bound,plot=1 or 0}
    //// Cutting maxAbs: {'CutAbs',maxAbs, action: warning = 0 or cutting = 1 ,plot=1 or 0}
    //// example: 
    //pretreat = list({'MA',3,1,1}, {'Deriv',1,1},{'MA',1,2,1},{'Cut',240,350,1});
    ////pretreat = list({'MA',3,1,1}, {'Loess',0.2,1,1},{'Cut',240,350,1},{'CutAbs',3.5,1,1});
    //
    //// Set of Analysis
    //// Lambert-Beer: {'LB'}
    //// Principal Component Analysis (PCA): {'PCA'} 
    //analysis = list({'LB'},{'PCA'});
    ////analysis = list({'PCA'})
    ////analysis = []
    //
    //// outlier analysis
    //outlier = 1; // 0 = no; 1 = yes.
    //
    //// Gravar Xval para os diferentes número de regressores
    //gravarXval = 1; //1 = grava, 0 = não grava
    //
    //
    //==========================================
    //   Pretreatment
    //==========================================
    lambda0=absor0(1,:); // first row is the wavenumber (or wavelength)
    absor0=absor0(2:$,:);
    absor = absor0;
    lambda=lambda0;
    x = x0;
    ifig = 1;

    [absor,lambda,x,ifig]=func_pretreatment(pretreat,absor,lambda,x,ifig)

    [absor,lambda,x,ifig,nclusters,idclusters]=func_analysis(analysis,absor,lambda,x,ifig)

    xtest=dadosteste{1}
    if ~isempty(xtest) then
        absortest=dadosteste{2}
        lambdatest = absortest(1,:)
        absortest = absortest(2:$,:)
        [absortest,lambdatest,xtest,ifig]=func_pretreatment(pretreat,absortest,lambdatest,xtest,ifig)
    end
    
    /////////////////////////////////////////////////////////////////////////////////////////////

    disp('-----------------------'); disp('  ');

    //-// Teste de Erro 1 - dimensoes das matrizes
    [ndx,nlx]=size(x);
    [nd,nl]=size(absor);
    if nd ~= ndx  // dados experimentais
        disp(['número de dados do arquivo de Concentração (=' string(ndx) ...
        ') não confere com o do arquivo de Absorbâncias (=' string(nd) ')'])
        disp('Programa Finalizado')
        return;
    end
    if nlx ~= nc
        disp(['número de componentes (=' string(nc) ') não confere com o arquivo dado (=' string(nlx) ')'])
        disp('Programa Finalizado')
        return;
    end



    //// Calculando concentrações máximas e mínimas
    xmax = max(x,'r');
    xmin = min(x,'r');
    // Calculando absorbâncias máximas e mínimas/ para cada coluna, ou seja,
    // para cada lambda
    Amax = max(absor,'r');
    Amin = min(absor,'r');



    // Normalizando as concentraçoes de 0 a 1
    x = x./repmat(xmax,nd,1)


    //// Separando dados para teste


    if ~isempty(xtest) then
        ntest = size(xtest,1)
        xtest =xtest./repmat(xmax,ntest,1)
    elseif frac_test >0 // só entra se não houver dados externos
        nd = size(x,1);
        //ind = [1:nd]';
        //lixo = sortrows([rand(nd,1) ind],1);
        [lixo,ind]=gsort(rand(nd,1));
        ind = ind(1:floor(nd*frac_test));
        xtest = x(ind,:);
        x(ind,:)=[];
        absortest = absor(ind,:);
        absor(ind,:)=[];
        ntest = size(xtest,1);
    end
    
    [nd,nl]=size(absor);    //nd = no. de dados experimentais, nl = no. compr. onda.

    // Model complexity evaluation
    select OptimModel{1}
    case 'Val' then
        frac_val = OptimModel{2}
        if frac_val <=0 then; disp('frac_val <=0'); abort; end
        [lixo,ind]=gsort(rand(nd,1));
        ind = ind(1:floor(nd*frac_val));
        Xval = x(ind,:);
        xaj = x;
        xaj(ind,:)=[];
        absorval = absor(ind,:);
        absoraj=absor;
        absoraj(ind,:)=[];
        nval = size(Xval,1);
        [ndaj,nl]=size(absoraj);
    case 'kfold' then
        kpart = OptimModel{2}
        if kpart == -1 then
            kpart = nd //LOOCV
            ndaj = nd-1 // data for fitting
        else
            [lixo,indval]=gsort(rand(nd,1)); // data shuffling
            x = x(indval,:) 
            absor = absor(indval,:)
            npart = ceil(nd/kpart) // the last partition may have less data then npart 
            ndaj = nd-npart // data for fitting

        end
        Xval = x
    case 'IC' then
        kpart = 1
        ndaj = nd
        select OptimModel{2}
        case 'AIC'
            rho = OptimModel{3}
        end
    else
        disp('Wrong OptimModel')
        abort
    end


    //// Evitando a singularida das matrizes
    if nl > ndaj //matrizes sigulares
        kmaxmax = ndaj;
        disp('O número de regressores possíveis do modelo está limitado pelo número de dados experimentais!');
    else
        kmaxmax = nl;
        disp('O número de regressores possíveis do modelo está limitado pelo número de comprimentos de onda!');
    end
    disp('Continuando . . .' )

    if kmax > kmaxmax
        kmax = kmaxmax;
    end




    //// Escolha lambdas SPA

    if Selecao == 2
        if optkini == 2 // otimiza lini
            for ind0 = 1:nl
                ls0(:,ind0) = spa_clean(absor,ind0,kmax)
            end
            conta=zeros(nl,1);
            for ind0 = 1:nl
                inds = find(ls0==ind0);
                conta(ind0) = sum(inds); 
            end
            [lixo,ind]=gsort(conta);
            cini = ind(1); // coluna escolhida inicialmente
            lini = lambda(cini);

        elseif optkini ==1
            cini = find(lambda == lini);
            if isempty(lixo)
                disp('---------')
                disp('O comprimento de onda escolhido para o início do SPA não está na região permitida!')
                disp('Comprimento de onda escolhido')
                disp(lini)
                disp('=> Finalizado')
                return
            end
        else
            lini = lambda(1);
            cini = 1
        end
        disp('lini')
        disp(lini)
    end




    //// Validação cruzada
    select OptimModel{1}
    case 'Val' then
        Xpval = zeros(Xval);
    case 'kfold' then
        Xpval = zeros(nd,nc);
    case 'IC' then
        IC = zeros(kmax,nc)
        Xpval = [];
    end

    SEcal=zeros(kmax,nc)
    lado = ceil(sqrt(nc));
    Xp = zeros(ndaj,nc);
    Xpp = zeros(nd,nc);
    Xpt = zeros(xtest);
    RMSEtest = zeros(kmax,nc);
    RMSECV = zeros(kmax,nc);
    RMSEcal = zeros(kmax,nc);


    formato0= '';
    formato= '';
    for i=1:size(cname,2)-1
        formato0 = formato0 + '%s \t';
        formato = formato + '%f \t';
    end
    formato0 = formato0 +' %s \n';
    formato = formato +' %f \n';


    for k=1:kmax //número máximo de regressores a ser testado
        disp(' k' ); disp(k);
        scf(ifig);
        for j=1:nc
            select OptimModel{1}
            case 'Val' then
                /////////////////
                select  Selecao
                case 1 then
                    [Xp(:,j),Xpval(:,j),Par_norm] = pls_model(absoraj,xaj(:,j),k,absorval)
                case 2 then
                    [Xp(:,j),Xpval(:,j),Par_norm] = spa_model(absoraj,xaj(:,j),k,cini,absorval)
                case 3 then
                    [Xp(:,j),Xpval(:,j),Par_norm] = pcr_model(absoraj,xaj(:,j),k,absorval)
                case 4 then
                    //lolimot parametros
                    data = [absoraj xaj(:,j)];sigma=0.33;nbpart=k;maximp=0;nbCut=2;vec=0;Log=1;modelIn=0;pLinear=0;
                    [modelOut, stat] = lolimot_learn(data, sigma, nbpart, maximp, nbCut, vec, Log)
                    for aux = [1:size(data,'r')]
                         Xp(aux,j)   = lolimot_estim(absoraj(aux,:), modelOut)
                    end
                    for aux=[1:size(absorval,'r')]
                         Xpval(aux,j)= lolimot_estim(absorval(aux,:), modelOut)
                    end
                else
                    disp('Erro - Escolha um algoritmo')
                    return
                end
                //////////////
                RMSECV(k,j)=sqrt(sum((Xpval(:,j)-Xval(:,j)).^2)/nval);
                RMSEcal(k,j)=sqrt(sum((Xp(:,j)-xaj(:,j)).^2)/ndaj);

            case 'kfold' then
                i1=0
                naj = 0
                
                for i=1:kpart
                    ////////////////
                    i0 = i1+1
                    i1 = min(i0+npart-1,nd)
                    naj = naj +i1-i0+1
                    absoraj = absor;
                    xaj=x;
                    absorval = absor(i0:i1,:);
                    absoraj(i0:i1,:) = []; //retirando linha de validação
                    xaj(i0:i1,:)=[]; // retirando a concentração
                    Xp = zeros(xaj); // pode variar
                    select  Selecao
                    case 1 then
                        [Xp(:,j),Xpval(i0:i1,j),Par_norm] = pls_model(absoraj,xaj(:,j),k,absorval)
                    case 2 then
                        [Xp(:,j),Xpval(i0:i1,j),Par_norm] = spa_model(absoraj,xaj(:,j),k,cini,absorval)
                    case 3 then
                        [Xp(:,j),Xpval(i0:i1,j),Par_norm] = pcr_model(absoraj,xaj(:,j),k,absorval)
                    case 4 then
                        data = [absoraj xaj(:,j)];sigma=0.33;nbpart=k;maximp=0;nbCut=2;vec=0;Log=1;modelIn=0;pLinear=0;
                        [modelOut, stat] = lolimot_learn(data, sigma, nbpart, maximp, nbCut, vec, Log)
                        for aux = [1:size(data,'r')]
                             Xp(aux,j)   = lolimot_estim(absoraj(aux,:), modelOut)
                        end
                        for aux=[1:size(absorval,'r')]
                             Xpval(i0-1+aux,j)= lolimot_estim(absorval(aux,:), modelOut)
                        end
                    else
                        disp('Erro - Escolha um algoritmo')
                        return
                    end
                    ///////////////
                    SEcal(k,j)= SEcal(k,j)+sum((Xp(:,j)-xaj(:,j)).^2);
                end
                // ajuste com todos para checar teste
                RMSECV(k,j)=sqrt(sum((Xpval(:,j)-Xval(:,j)).^2)/nd);
                RMSEcal(k,j)=sqrt(SEcal(k,j)/naj);
                //        RMSEcal(k,j)=sqrt(SEcal(k,j)/(nd*(nd-1-k))); // desvio padrão de ajuste
            case 'IC' then
                /////////////////////
                absoraj = absor;
                xaj=x;
                select  Selecao
                case 1 then 
                    [Xp(:,j),lixo,Par_norm] = pls_model(absoraj,xaj(:,j),k,[])
                case 2 then
                    [Xp(:,j),lixo,Par_norm] = spa_model(absoraj,xaj(:,j),k,cini,[])
                case 3 then
                    [Xp(:,j),lixo,Par_norm] = pcr_model(absoraj,xaj(:,j),k,[])
                case 4 then
                        data = [absoraj xaj(:,j)];sigma=0.33;nbpart=k;maximp=0;nbCut=2;vec=0;Log=1;modelIn=0;pLinear=0;
                        [modelOut, stat] = lolimot_learn(data, sigma, nbpart, maximp, nbCut, vec, Log)
                        for aux = [1:size(data,'r')]
                             Xp(aux,j)   = lolimot_estim(absoraj(aux,:), modelOut)
                        end
                else
                    disp('Erro - Escolha um algoritmo')
                    return
                end
                //////////////////
                RMSEcal(k,j)=sqrt(sum((Xp(:,j)-xaj(:,j)).^2)/ndaj);
                select OptimModel{2}
                case 'AIC'
                    IC(k,j) = ndaj*log(RMSEcal(k,j)^2/0.05^2) + rho*k
                end
            end


            if ~isempty(xtest)
                switch  Selecao
                /////////////////
                case 1
                    [Xpp(:,j),Xpt(:,j),Par_norm] = pls_model(absor,x(:,j),k,absortest)
                case 2
                    [Xpp(:,j),Xpt(:,j),Par_norm] = spa_model(absor,x(:,j),k,cini,absortest)
                case 3
                    [Xpp(:,j),Xpt(:,j),Par_norm] = pcr_model(absor,x(:,j),k,absortest)
                case 4 then
                    data = [absor,x(:,j)];sigma=0.33;nbpart=k;maximp=0;nbCut=2;vec=0;Log=1;modelIn=0;pLinear=0;
                    [modelOut, stat] = lolimot_learn(data, sigma, nbpart, maximp, nbCut, vec, Log)
                    for aux=[1:size(absor,'r')]
                        Xpp(aux,j)   = lolimot_estim(absor(aux,:), modelOut)
                    end
                    for aux=[1:size(absortest,'r')]
                        Xpt(aux,j)= lolimot_estim(absortest(aux,:), modelOut)
                    end
                otherwise
                    disp('Erro - Escolha um algoritmo')
                    return
                end
                /////////////////////////
                RMSEtest(k,j) = sqrt(sum((Xpt(:,j)-xtest(:,j)).^2)/ntest);
                RMSEajtest(k,j) = sqrt(sum((Xpp(:,j)-x(:,j)).^2)/nd);
                
                if OptimModel{1} == 'IC'
                    Xval = x
                    Xpval= Xp
                end  
                set(gca(),"auto_clear","on")
                subplot(lado,lado,j); plot(Xval(:,j),Xpval(:,j),'o',xtest(:,j),Xpt(:,j),'x',[0 1], [0 1], '-');
                legend( 'val', 'test',-1 );
            else
                if OptimModel{1} == 'IC'
                    Xval = x
                    Xpval= Xp
                end
                set(gca(),"auto_clear","on")
                subplot(lado,lado,j); plot(Xval(:,j),Xpval(:,j),'o',[0 1], [0 1], '-');
                legend('val',-1);
            end

            xlabel = ('Cref');ylabel = ('Cpred');
            title = (['no. regressores = ',  string(k), ' componente = ', cname(j)]);
            xtitle(title,xlabel,ylabel);
        end

        // outlier analysis
        if outlier == 1 then
            Erroval = Xpval-Xval
            Std_loo = zeros(Erroval)
            for i=1:nd
                err_loo = Erroval;
                err_loo(i,:) = [];
                Std_loo(i,:) = sqrt(diag(err_loo'*err_loo)/(nd-1))';
            end
            tcalc = Erroval./Std_loo
            [Pstat,Qstat]=cdft("PQ",abs(tcalc),(nd-1)*ones(tcalc));
            [indr,indc]=find(Pstat>.995); // nível de confiança de 99%;
            if ~isempty(indr) then
                indr=unique(indr);
                disp('outliers:');
                disp(['line', 'p-value < 0.005']);
                disp([indr' Qstat(indr,:) ])
            end
        end
        if gravarXval == 1 then
            // Gravando Xval para os diferentes modelos
            
            Xvalconc(indval,:) = Xpval.*repmat(xmax,nd,1)
            [fd, err] = mopen('Xvalconc_'+string(k)+'.txt' , 'w')
            mfprintf(fd,formato0,cname);
            mfprintf(fd,formato,Xvalconc);
            mclose(fd)
        end
    end
    ifig=ifig+1;

    //// RESULTADOS

    if  OptimModel{1} == 'IC'
        RMSECV = IC
    end
    for i=1:nc
        // em unidades de concentração =unidades do arq de dados
        RMSECVconc(:,i)=RMSECV(:,i)*xmax(i);  
        RMSEcalconc(:,i)=RMSEcal(:,i)*xmax(i);
        if ~isempty(xtest) then
            RMSEtestconc(:,i)=RMSEtest(:,i)*xmax(i);
            scf(ifig); plot([1:kmax]',RMSECV(:,i),[1:kmax]',RMSEcal(:,i),[1:kmax]',RMSEtest(:,i),'--',[1:kmax]',RMSEajtest(:,i),'--');
            legend('RMSECV','RMSEcal','RMSEtest','RMSEajtest')
        else
            scf(ifig); plot([1:kmax]',RMSECV(:,i),[1:kmax]',RMSEcal(:,i));
            legend('RMSECV','RMSEcal')
        end
        title = cname(i); ifig=ifig+1;
        xlabel = ('Latent Variables');ylabel = ('RMSECV');
        xtitle(title,xlabel,ylabel);
    end


    //format short e
    //format short
    cnametxt = ' ';
    for i=1:size(cname,2)//length(cname)
        //cnametxt = cnametxt + '   ' + cname{i};
        cnametxt = cnametxt + '   ' + cname(i);
    end
    disp('')
    disp('******************************************')
    disp('')
    disp(['k', cnametxt])
    disp(' RMSECVn')
    disp([[1:kmax]' RMSECV])
    disp(' RMSECVconc')
    disp([[1:kmax]' RMSECVconc])
    disp('')
    disp(['k', cnametxt])
    disp(' RMSEcaln')
    disp([[1:kmax]' RMSEcal])
    disp(' RMSEcalconc')
    disp([[1:kmax]' RMSEcalconc])

    if ~isempty(xtest) then
        disp(' RMSEtestn')
        disp([[1:kmax]' RMSEtest] )
        disp(' RMSEtestconc')
        disp([[1:kmax]' RMSEtestconc] )
    end
    disp('')
    disp('')
    disp('******************************************')


    formato0= 'k \t';
    formato= '%d \t';
    for i=1:size(cname,2)-1
        formato0 = formato0 + '%s \t';
        formato = formato + '%f \t';
    end
    formato0 = formato0 +' %s \n';
    formato = formato +' %f \n';


    [fd, err] = mopen('Erro_cv.txt' , 'w')
    mfprintf(fd,formato0,cname);
    mfprintf(fd,formato,[[1:kmax]' RMSECVconc]);
    mclose(fd)
    if ~isempty(xtest) then
        [fd, err] = mopen('Erro_test.txt' , 'w')
        mfprintf(fd,formato0,cname);
        mfprintf(fd,formato,[[1:kmax]' RMSEtestconc]);
        mclose(fd)
    end

    if Selecao == 2 then
        disp('comprimentos de onda escolhidos')
        disp(lambda(Par_norm.ind));
    end
    flag = 0; // saída
endfunction

