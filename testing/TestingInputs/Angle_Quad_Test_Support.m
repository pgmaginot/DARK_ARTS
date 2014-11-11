n_dir = 8;

[mu,w] = GLNodeWt(8);

for d=1:1:n_dir
   fprintf(' expected_mu[%i] = %12.8f ; expected_w[%i] = %12.8f ; \n',d-1,mu(d,1) , d-1 , w(d,1) ) 
end