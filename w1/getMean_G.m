function rls = getMean_G(price)

 rls = prod((price(2:end)./price(1:end-1))) ^ (1/(length(price))) - 1;
 

