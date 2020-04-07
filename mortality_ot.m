%mortality_exchange
% Based on the mortality of DM I, the mortality of the same amount of DM II can be calcualted.
% For Anopheles, it was assume that one unit of DM II was equal to twelve units of DM I.
% the slope of the probit-dose curve is 2, based on the probit-dose curve reported by previous works.
function mortality_out=mortality_ot(IR)
mortality=IR;
if IR==0
    mortality=0.01;
end
M=mortality;
slope=2;
Probit_50=0;
LD50=61;
inter=Probit_50-slope*log10(LD50);
Probit=2^0.5*erfinv(2*M-1);
dose=10^((Probit-inter)/slope);
Insecticide_enhancment=12; 
doseout=Insecticide_enhancment*dose;
Probitout=inter+slope*log10(doseout);
mortality_out=(1+erf(Probitout/2^0.5))/2;
%%
if IR==1
    mortality_out=1;
end
end
