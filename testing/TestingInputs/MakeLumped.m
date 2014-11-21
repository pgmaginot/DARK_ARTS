function M = MakeLumped(M,n)

for i=1:1:n
    for j=1:1:n
        if(i==j)
            continue
        end
        M(i,i) = M(i,i) + M(i,j);
        M(i,j) = 0;
    end
end

return
end