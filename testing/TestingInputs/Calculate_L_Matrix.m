function [L,f] = Calculate_L_Matrix(is_positive , n_dfem , n_quad , q_w, dfem , d_dfem , dfem_edge) 

    L = zeros(n_dfem);
    f = zeros(n_dfem,1);

    for i=1:1:n_dfem
        if(is_positive)
            f(i,1) = dfem_edge(1,i);
        else
            f(i,1) = -dfem_edge(2,i);
        end

        for j=1:1:n_dfem
           % the only quadrature evaluations guaranteeed to have endpoint
           % evaluations are dat.b_f and dat.b_leak
            if(is_positive)
                L(i,j) = dfem_edge(2,j)*dfem_edge(2,i);
            else
                L(i,j) = -dfem_edge(1,j)*dfem_edge(1,i);
            end

            for p=1:1:n_quad
                if(is_positive)
                    L(i,j) = L(i,j) - q_w(p)*dfem(p,j)*d_dfem(p,i);
                else
                    L(i,j) = L(i,j) - q_w(p)*dfem(p,j)*d_dfem(p,i);
                end
            end           
        end        
    end
       
    return
end