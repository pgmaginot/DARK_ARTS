function Calculate_Matrices()

clc

fprintf('Exact_Integration_Test\n\n');
n_dfem = 3;
[x,w] = NewtonCotes(n_dfem);

n_q = n_dfem + 3;
[int_q,int_w] = GLNodeWt(n_q);

[dfem_evals, d_dfem_evals] = feshpln(int_q , x , n_dfem-1);
[dfem_edge,d_dfem_edge] = feshpln([-1;1],x,n_dfem-1);

M = GenerateMassMatrix( n_dfem , n_q , dfem_evals , int_w , false);
[Lp, fp ] = Calculate_L_Matrix(true , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
[Ln, fn ] = Calculate_L_Matrix(false , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
print_matrices(M,Lp,Ln , fp, fn,n_dfem);

clear all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nGauss_Self_Lumping_Test\n\n');
n_dfem = 5;
[x,w] = GLNodeWt(n_dfem);

n_q = n_dfem;
[int_q,int_w] = GLNodeWt(n_q);

[dfem_evals, d_dfem_evals] = feshpln(int_q , x , n_dfem-1);
[dfem_edge,d_dfem_edge] = feshpln([-1;1],x,n_dfem-1);
M = GenerateMassMatrix( n_dfem , n_q , dfem_evals , int_w , false);
[Lp, fp ] = Calculate_L_Matrix(true , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
[Ln, fn ] = Calculate_L_Matrix(false , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
print_matrices(M,Lp,Ln , fp, fn,n_dfem);

clear all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nLobatto_Self_Lumping_Test\n\n');
n_dfem = 4;
[x,w] = lglnodes(n_dfem);

n_q = n_dfem;
[int_q,int_w] = lglnodes(n_q);

[dfem_evals, d_dfem_evals] = feshpln(int_q , x , n_dfem-1);
[dfem_edge,d_dfem_edge] = feshpln([-1;1],x,n_dfem-1);
M = GenerateMassMatrix( n_dfem , n_q , dfem_evals , int_w , false);
[Lp, fp ] = Calculate_L_Matrix(true , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
[Ln, fn ] = Calculate_L_Matrix(false , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
print_matrices(M,Lp,Ln , fp, fn,n_dfem);

clear all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nTraditional_Lumping_Test\n\n');
n_dfem = 4;
[x,w] = NewtonCotes(n_dfem);

n_q = n_dfem + 2;
[int_q,int_w] = GLNodeWt(n_q);

[dfem_evals, d_dfem_evals] = feshpln(int_q , x , n_dfem-1);
[dfem_edge,d_dfem_edge] = feshpln([-1;1],x,n_dfem-1);
M = GenerateMassMatrix( n_dfem , n_q , dfem_evals , int_w , true);
[Lp, fp ] = Calculate_L_Matrix(true , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
[Ln, fn ] = Calculate_L_Matrix(false , n_dfem , n_q , int_w , dfem_evals , d_dfem_evals , dfem_edge) ;
print_matrices(M,Lp,Ln , fp, fn,n_dfem);



return
end