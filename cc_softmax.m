function y = cc_softmax(x,beta)

y = exp( beta*x) / sum( exp(beta*x)); 


end
