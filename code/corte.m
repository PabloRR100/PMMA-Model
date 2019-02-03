function f = corte(i,n)

    m = max(1,round((n-1)*rand));
    f(1,1) = m;
    
    j = round((i*(n-m)+m)/n);
    
        f(1,2)= i-j;
        f(2,1)= n-m;
        f(2,2)= j-1; % -1 El que se rompe

end

