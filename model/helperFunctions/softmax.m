function pAction = softmax(reward, rho, beta, bias)

pStay = 1/(1 + exp(bias-beta*(reward-rho))); 
pLeave = 1-pStay;
pAction = [pLeave pStay];

end



