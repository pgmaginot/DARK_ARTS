function print_matrices(M,Lp,Ln , fp, fn,n_dfem)

for i=1:1:n_dfem
    for j=1:1:n_dfem
        fprintf('expected_dimensionless_mass[%i][%i] = %10.8f ; \n',i-1,j-1,M(i,j) );
    end    
end

for i=1:1:n_dfem
    for j=1:1:n_dfem
        fprintf('expected_L_mu_positive[%i][%i] = %10.8f ; \n',i-1,j-1,Lp(i,j) );
    end    
end

for i=1:1:n_dfem
    for j=1:1:n_dfem
        fprintf('expected_L_mu_negative[%i][%i] = %10.8f ; \n',i-1,j-1,Ln(i,j) );
    end    
end

for i=1:1:n_dfem
    fprintf('expected_f_mu_positive[%i] = %10.8f ; \n',i-1,fp(i,1) ); 
end

for i=1:1:n_dfem
    fprintf('expected_f_mu_negative[%i] = %10.8f ; \n',i-1,fn(i,1) ); 
end


return
end