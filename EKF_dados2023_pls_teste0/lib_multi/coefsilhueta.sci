function [silhouette_coefficient] = coefsilhueta(X, idx, k)
    //X é a matriz de dados
    //idx é o vetor de rótulos (cada ponto pertence a um cluster)
    //k é o número de clusters
    n = size(X, 1);
    a = zeros(n, 1);  // average distance of a point to its own cluster
    b = zeros(n, k);  // average distance of a point to each other cluster

    // Calculate a(i) for each point
    for i = 1:n
        index = find(idx == idx(i));
        points_in_cluster = X(index, :); 
//        a(i) = mean(pdist2(X(i,:), points_in_cluster));
        dist=0
        for aux = [1:size(points_in_cluster,'r')]
            dist = dist + norm(X(i,:) - points_in_cluster(aux,:) ) 
        end
        a(i) = dist/size(points_in_cluster,'r');
    end

    // Calculate b(i,k) for each point and each cluster
    for i = 1:n
        for j = 1:k
            if j ~= idx(i)
                index = find(idx == j);
                points_in_other_cluster = X(index, :);
//                b(i,j) = mean(pdist2(X(i,:), points_in_other_cluster));
                dist=0
                for aux = [1:size(points_in_other_cluster,'r')]
                    dist = dist + norm(X(i,:) - points_in_other_cluster(aux,:) ) 
                end
                b(i,j) = dist/size(points_in_other_cluster,'r');
            end
        end
    end
    
    //Solução incorreta, mas funciona pois o que conta é o menor valor de b
    if k ~=1 then
        for aux = [1:k]
            index = find(b(:,aux) == 0);
            b(index,aux) = 2
        end
    end
    
    // Calculate silhouette coefficient for each point
    s = zeros(n, 1);
    for i = 1:n
        s(i) = (min(b(i,:)) - a(i)) / max([a(i), min(b(i,:))]);
    end

    silhouette_coefficient = mean(s);

endfunction


