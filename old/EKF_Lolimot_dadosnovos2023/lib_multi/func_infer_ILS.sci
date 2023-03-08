function [Xp,Xinf,Par_norm,absor,x,Xinf_conc,RMSECV] = func_infer_ILS(Selecao,optkini,lini,kinf,nc,cname,unid,x0,absor0,xinf0,absorinf0,pretreat,pretreatinf,analysis,analysisinf,leverage)

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
    flag = -1;
    ////==========================================================================================================////
    ////                                   INICIO DE INSERÇAO DE DADOS                                                        ////
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
    //
    ////// Método
    //Selecao = 3; // 1 = PLS; 2 = SPA; 3=PCR; 4=lolimot
    
    //==========================================
    //  Initializing Pretreatment
    //==========================================
    lambda0=absor0(1,:); // first row is the wavelength (or wavenumber)
    absor0=absor0(2:$,:);
    absor = absor0;
    lambda=lambda0;
    x = x0;
    ifig = 1;

    [absor,lambda,x,ifig]=func_pretreatment(pretreat,absor,lambda,x,ifig)

    [absor,lambda,x,ifig]=func_analysis(analysis,absor,lambda,x,ifig)


    lambdainf = absorinf0(1,:);
    absorinf=absorinf0(2:$,:);
    xinf = xinf0;


    [absorinf,lambdainf,xinf,ifig]=func_pretreatment(pretreatinf,absorinf,lambdainf,xinf,ifig)

    [absorinf,lambdainf,xinf,ifig]=func_analysis(analysisinf,absorinf,lambdainf,xinf,ifig)
 

    [nd,nl]=size(absor);    //nd = no. de dados experimentais, nl = no. compr. onda.
    [ndx,ncx]=size(x);
    [ndinf,nlinf]=size(absorinf);    //nd = no. de dados experimentais, nl = no. compr. onda.

    if ~isempty(xinf)
        [ndxinf,ncxinf]=size(xinf);
    end

    if nd ~= ndx  // dados experimentais
        disp(['número de dados do arquivo de Concentração (=' string(ndx) ...
        ') não confere com o do arquivo de Absorbâncias (=' string(nd) ')'])
        disp('Programa Finalizado')
        return;
    end
    if ncx ~= nc
        disp(['número de componentes (=' string(nc) ') não confere com o arquivo dado (=' string(ncx) ')'])
        disp('Programa Finalizado')
        return;
    end
    if ~isempty(xinf)
        if ndinf ~= ndxinf  // dados experimentais
            disp(['número de dados do arquivo de Concentração inf(=' string(ndxinf) ...
            ') não confere com o do arquivo de Absorbâncias inf (=' string(ndinf) ')'])
            disp('Programa Finalizado')
            return;
        end
    end
    if ~isempty(xinf)
        if ncxinf ~= nc
            disp(['número de componentes inferidos (=' string(nc) ') não confere com o arquivo dado (=' string(ncxinf) ')'])
            disp('Programa Finalizado')
            return;
        end
    end

    if nl ~= nlinf
        disp(['número de lambdas (=' string(nl) ') não confere com o arquivo dado (=' string(nlinf) ')'])
        disp('Programa Finalizado')
        return;
    end


    disp('Pausa');
    disp('escreva ''resume'' ou ''abort''');
    pause; disp('Continuando...')



    //// Calculando concentrações máximas e mínimas
    xmax = max(x,'r');
    xmin = min(x,'r');
    // Calculando absorbâncias máximas e mínimas/ para cada coluna, ou seja,
    // para cada lambda
    Amax = max(absor,'r');
    Amin = min(absor,'r');



    // Normalizando as concentraçoes de 0 a 1
    x = x./repmat(xmax,nd,1)


    //// Evitando a singularida das matrizes
    if nl > nd //matrizes sigulares
        kmaxmax = nd;
        disp('O número de regressores possíveis do modelo está limitado pelo número de dados experimentais!');

    else
        kmaxmax = nl;
        disp('O número de regressores possíveis do modelo está limitado pelo número de comprimentos de onda!');

    end
    disp('Continuando . . .' )

    if kinf > kmaxmax
        disp('erro -- kinf > máximo');
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

    //// Inferência
    for j=1:nc
        k = kinf(j);
        absoraj = absor;
        xaj=x;
        select  Selecao
        case 1 then
            [Xp(:,j),Xinf(:,j),Par_norm] = pls_model(absoraj,xaj(:,j),k,absorinf)
        case 2 then
            [Xp(:,j),Xinf(:,j),Par_norm] = spa_model(absoraj,xaj(:,j),k,cini,absorinf)
        case 3 then
            [Xp(:,j),Xinf(:,j),Par_norm] = pcr_model(absoraj,xaj(:,j),k,absorinf)
        case 4 then
            data = [absoraj,xaj(:,j)];sigma=0.33;nbpart=k;maximp=0;nbCut=2;vec=0;Log=1;modelIn=0;pLinear=0;
            [modelOut, stat] = lolimot_learn(data, sigma, nbpart, maximp, nbCut, vec, Log)
            for aux=[1:size(absoraj,'r')]
                Xp(aux,j)   = lolimot_estim(absoraj(aux,:), modelOut)
            end
            for aux=[1:size(absortest,'r')]
                Xinf(i0-1+aux,j)= lolimot_estim(absoraj(aux,:), modelOut)
            end    
        else
            disp('Erro - Escolha um algoritmo')
            return
        end
        
        // Leverage analysis
        if leverage == 1 then
            disp('Evaluating outlier: leverage - '+cname{j})
            select Selecao
            case 2 then
                T = absor(Par_norm.ind)
            else
                T = Par_norm.T;
            end
            lT = size(T,1)
            meanT = mean(T,'r') ;
            ET = T - repmat(meanT,lT,1);
            hii = 1/lT + diag(ET*inv(T'*T)*ET');
            hcrit = 3*k/lT;
            ind_lev = find(hii>hcrit);
            if ~isempty(ind_lev) then
                disp('leverage - hii > hcrit para as linhas')
                disp(ind_lev')
            end
        end
    end

    //// RESULTADOS
    ninf = size(Xinf,1)
    Xinf_conc = Xinf.*repmat(xmax,ninf,1)
    if ~isempty(xinf)  then
        for i=1:nc
            scf(10+i); plot(xinf(:,i),Xinf_conc(:,i),'o',[min(xinf(:,i)),max(xinf(:,i))],[min(xinf(:,i)),max(xinf(:,i))],'-')    
            xtitle(cname(i)+' ('+ unid + ')','Conc exp', 'Conc inferida')
        end
    end
    disp('Concentrações Preditas')
    texto = [] ;
    for i=1:size(cname,2)
        texto = [texto cname(i) ] ;
    end
    disp(texto)
    disp(Xinf_conc)

    /// Gravando resultado em arquivo de dados

    formato= '%f \t';
    for i=2:size(Xinf_conc,2)-1
        formato = formato + '%f \t';
    end
    formato = formato +' %f \n';
    [fd, err] = mopen('Xinf.txt' , 'w')
    mfprintf(fd,formato,Xinf_conc);
    mclose(fd)

    flag = 0;
   
endfunction
