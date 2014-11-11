function Evaluating_Legendre_Polynomials(x_eval)

p_ord = 4;

poly_eval = zeros(5,1);


for i=0:1:p_ord
    temp = legendre( i , x_eval); 
    poly_eval(i+1,1) = temp(1,1);
    fprintf('%12.8f\n',poly_eval(i+1,1) )
end

end