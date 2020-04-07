function p = DMParameters()

p.muE = 0.2;

% Gompertz for mortality curve
p.muA = [0.0062    0.0074    0.0089    0.0107    0.0129    0.0155    0.0186    0.0223    0.0268    0.0321 ...
    0.0384    0.0461    0.0551    0.0659    0.0787    0.0939    0.1118    0.1328    0.1575    0.1862 ... 
    0.2195    0.2576    0.3011    0.3500    0.4043    0.4636    0.5271    0.5937    0.6614    0.7280 ...
    0.7911    0.8478    0.8960    0.9343    0.9621    0.9805    1];

p.dmin = 0.035;
p.dmax = 0.4;
p.K = 5e5;
p.c1 = 0.04939737836602;
p.c2 = 3.3978873097;
p.c3 = 2.0;

p.timesteps = 365;%simulation for 3650 days, which are 10 years.

p.delay = 12;% mosquitoes required at least 12 days following an infectious blood meal to become infectious to humans.

p.b = 0.55;
p.recovery = 1-exp(-1/25); 

p.k1 = 0.02;
p.k2 = 0.2;

p.batch = 60;
p.batchExp = 60;
p.female = 0.5;
p.feedfirst = 0.5;
p.matefirst = 1 - p.feedfirst;

p.insect = 0; % reduction in survival rate for insecticide %the mortality of insects after contact insecticides.
p.hs = 0; % reduction in susceptibility to compound

