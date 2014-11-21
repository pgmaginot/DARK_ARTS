function M_quad = GenerateMassMatrix(n_dfem,n_qp,b,w,lump)


M_quad = zeros(n_dfem);
for i=1:1:n_dfem
    for j=1:1:n_dfem
        for p=1:1:n_qp
            M_quad(i,j) = M_quad(i,j) + w(p)*b(p,i)*b(p,j);
        end
    end
end
        
if(lump)
    M_quad = MakeLumped(M_quad,n_dfem);
end

return
end