  function [ sum_lik  ] = MLE_t5_tgarch( alfa_plus_hat,  alfa_minus_hat  ,betahat, sigma, y,h )



 T=2000;
 likelihoods= zeros(1,1999); 

 
 
 
for  j=2: T;
     h(j)= 1+ betahat .* h(j-1) +  alfa_plus_hat .* ((y(j-1)).^2) .* ((y(j-1))>0) + alfa_minus_hat .* ((y(j-1)).^2) .* ((y(j-1))<0)   ;
likelihoods(j) =  log(gamma(3)) - log(gamma(2.5))  - 0.5 .*log(pi) -0.5 .* log(5) - ...
                   3 .*    log(        1+ ( (     ( (y(j)^2) .* (5/3) )  ./ ( h(j) .* (sigma.^2)     )         )/5)     )    - ...
                   log(sigma ) - 0.5.* log(h(j)) + 0.5 .* log((5/3))   ; 

end
sum_lik= -sum(likelihoods) ;
 
  end
