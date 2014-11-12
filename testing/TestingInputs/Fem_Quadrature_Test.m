function Fem_Quadrature_Test()

n_dfem = 4;
n_xs = 5;
n_int = 8;

[dfem_pts , dfem_wts] = lglnodes( n_dfem );
[xs_pts , xs_wts] = NewtonCotes( n_xs );
[int_pts, int_wts] = GLNodeWt( n_int );

for i=1:1:n_dfem
   fprintf('dfem_pts[%i] =  %12.8f ;\n',i-1,dfem_pts(i,1)) 
end
Fem_Q
for i=1:1:n_xs
   fprintf('xs_pts[%i] =  %12.8f ;\n',i-1,xs_pts(i,1)) 
end

for i=1:1:n_int
   fprintf('integration_pts[%i] =  %12.8f ; \n',i-1,int_pts(i,1)) 
end

[dfem_at_int , d_dfem_at_int ] = feshpln(int_pts , dfem_pts , n_dfem - 1);
[xs_at_int , d_xs_at_int ] = feshpln(int_pts , xs_pts , n_xs - 1);


% we are storing all quadrature points for every basis function
cnt = 0;
for b=1:1:n_dfem
    for p=1:1:n_int
        fprintf('dfem_at_integration_pts[%i] = %12.8f ;\n',cnt , dfem_at_int(p,b))
    cnt = cnt + 1;
    end
end

cnt = 0;
for b=1:1:n_dfem
    for p=1:1:n_int
        fprintf('dfem_deriv_at_integration_pts[%i] = %12.8f ;\n',cnt , d_dfem_at_int(p,b))
    cnt = cnt + 1;
    end
end

cnt = 0;
for b=1:1:n_xs
    for p=1:1:n_int
        fprintf('xs_at_integration_pts[%i] = %12.8f ;\n',cnt , xs_at_int(p,b))
    cnt = cnt + 1;
    end
end



end