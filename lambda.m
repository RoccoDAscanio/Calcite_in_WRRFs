
function lambda_calc = lambda(S,IR,IZ,model_flag)

I = 0.019*S;
A = 0.509;

if strcmp(model_flag,'Davies')

    frac = -A*(IZ^2)*(((sqrt(I))/(1+sqrt(I)))-(0.3*I));

elseif strcmp(model_flag,'EDH')

    B = 0.328;
    frac = (-A*(IZ^2)*sqrt(I))/(1+(B*IR*sqrt(I)));

elseif strcmp(model_flag,'Ideal')

    frac = 0;

else

    error('Ion Activity Coefficient Model Flag Incorrectly Specified');

end

lambda_calc = 10^(frac);

end
