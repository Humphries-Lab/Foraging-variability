function PAction = epsilon_softmaxConstrain(valueStay, valueLeave, Beta, Epsilon)

numerator = exp([valueLeave, valueStay] * Beta); % softmax function
numerator(numerator==inf) = realmax/2; % replace any infs with just a large number
numerator(numerator==0) = eps(0); % replace any 0's with just a tiny number to prevent NaN emerging during division
PAction = numerator ./ sum(numerator) * (1-2*Epsilon) + Epsilon;

PAction(PAction == 0) = eps(0); % if the selected action has PAction = 0, then log(0) becomes infinity when calculating log likelihood (and breaks fmincon). Replace with a small number
end