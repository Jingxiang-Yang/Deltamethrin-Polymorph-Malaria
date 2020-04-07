%%
% On all feeding days mosquitoes will transmit disease. On the first resting day (R1) mosquitoes will contain insecticide but will not
% transmit disease.
%%
function [IHout,cov_vec] = DM(p)

cov_vec = p.cov_vec;
IHout = NaN(length(cov_vec),1);

for j=1:length(cov_vec)
    coverage = cov_vec(j);
    %% Pre-allocation
    LaidEggs = NaN(p.timesteps,1);
    E = NaN(p.timesteps,3);
    muL = NaN(p.timesteps,1);
    Ltot = NaN(p.timesteps,1);
    L = NaN(p.timesteps,10);
    R = NaN(p.timesteps,1);
    MateFirst = NaN(p.timesteps,36);
    FeedFirst = NaN(p.timesteps,36);
    MateFirstExp = zeros(p.timesteps,36);
    FeedFirstExp = zeros(p.timesteps,36);
    
    betaH = NaN(p.timesteps,1);
    betaM = NaN(p.timesteps,1);
    IH = NaN(p.timesteps,1);
    feeders = NaN(p.timesteps,1);
    
    %% Initial Conditions
    E(1,:) = 1;
    L(1,:) = 1;
    R(1,:) = 1;
    MateFirst(1,:) = 100;
    FeedFirst(1,:) = 100;
    IH(1) = .047;
    feeders(1) = 0;
    
    %% Running the Simulation
    
    for t = 1:p.timesteps-1
        %% Eggs
        LaidEggs(t) = p.batch * (sum(MateFirst(t,5:4:36))+sum(FeedFirst(t,5:4:36)))+ ...
            p.batchExp * (sum(MateFirstExp(t,5:4:36))+sum(FeedFirstExp(t,5:4:36)));
        
        E(t+1,1) = LaidEggs(t);
        E(t+1,2:end) = (1-p.muE)*E(t,1:end-1);
        
        %% Larvae/Pupae
        Ltot(t) = sum(L(t,:));
        muL(t) = p.dmin - p.c1 + (p.dmax-(p.dmin-p.c1))/(1+exp(-p.c2*Ltot(t)/p.K+p.c3));

        L(t+1,1) = (1-p.muE)*E(t,end);
        L(t+1,2:end) = (1-muL(t))*L(t,1:end-1);
        
        %% Adult stages
        R(t+1) = p.female*(1-muL(t))*L(t,end);
        
        MateFirst(t+1,1) = p.matefirst * (1-p.muA(1)) * R(t);
        MateFirst(t+1,2:end) = (1-p.muA(2:end-1)) .* MateFirst(t,1:end-1);
        MateFirstExp(t+1,2:end) = (1-p.muA(2:end-1)) .* MateFirstExp(t,1:end-1);
        % assumuption: For IRS, mosquitoes contact insecticide only on the indoor resting day, which is R1 day.
        MateFirstExp(t+1,3:4:36) = MateFirstExp(t+1,3:4:36) + (1-p.insect) * coverage * MateFirst(t+1,3:4:36);
        MateFirst(t+1,3:4:36) = MateFirst(t+1,3:4:36) - coverage * MateFirst(t+1,3:4:36);
        
        FeedFirst(t+1,1) = p.feedfirst * (1-p.muA(1)) * R(t);
        FeedFirst(t+1,2:end) = (1-p.muA(2:end-1)) .* FeedFirst(t,1:end-1);
        FeedFirstExp(t+1,2:end) = (1-p.muA(2:end-1)) .* FeedFirstExp(t,1:end-1);
        FeedFirstExp(t+1,3:4:36) = FeedFirstExp(t+1,3:4:36) + (1-p.insect)  *coverage * FeedFirst(t+1,3:4:36);
        FeedFirst(t+1,3:4:36) = FeedFirst(t+1,3:4:36) - coverage * FeedFirst(t+1,3:4:36);
        
        %% calcualte the number of infectious feeder
        %every feed behavior lead to possible infection.
        if t>p.delay
            feeders(t) = betaM(t-p.delay) *...
                (sum(MateFirst(t,2+p.delay:4:36)) +...
                sum(FeedFirst(t,[1+p.delay 6+p.delay:4:36]))) +...
                betaM(t-p.delay)*  ...
                (sum(MateFirstExp(t,2+p.delay:4:36)) +...
                sum(FeedFirstExp(t,[1+p.delay 6+p.delay:4:36])));
            betaH(t) = 1 - (1-p.b)^(p.bites*feeders(t));
        else
            betaH(t) = .1;
        end
        %% number of infected humans
        IH(t+1) = IH(t)+betaH(t)*(1-IH(t))-p.recovery*IH(t); % infected humans
        betaM(t) = p.k2 * (1-p.k1/(p.k2*IH(t)+p.k1)); % risk of becoming infeected for a mosquito
        
    end % end time loop
    
    IHout(j) = IH(end);
    
end % end coverage loop
