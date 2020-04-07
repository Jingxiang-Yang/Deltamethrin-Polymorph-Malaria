%% FeedersInfectedHuman_NoDM_DMI_DMII_45prevalence.m
%% Required functions:
% DM.m
% DMParameters.m
% mortality_ot.m

%%
%%45% PREVALENCE
% On all feeding days mosquitoes will transmit disease. Only on the first 
% resting day (R1) mosquitoes will contain insecticide but will not
% transmit disease but will not transmit disease.

%%
clc;clear;
%%
p = DMParameters();
coverage = 0.75;% the selected coverage
mortalityDMI=0.3;% corrsponding to 70% resistance.
p.bites = 55/p.K;
disp(p.bites)

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
IH(1,:) = .047;
feeders(1,:) = 0;
%% Running the Simulation
t_0=1;
t_DMI=p.timesteps*3*0.1;
t_DMII=p.timesteps*6*0.1;
t_end=p.timesteps*9*0.1;
%%
for t = 1:t_end
    if t==t_0
        mortality = 0; % no insecticide
    end
    if t==t_DMI
        mortality = mortalityDMI; % DM I
    end
    if t==t_DMII
        mortality=mortality_ot(mortality);% DM II
    end
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
    MateFirstExp(t+1,3:4:36) = MateFirstExp(t+1,3:4:36) + (1-mortality) * coverage * MateFirst(t+1,3:4:36);
    MateFirst(t+1,3:4:36) = MateFirst(t+1,3:4:36) - coverage * MateFirst(t+1,3:4:36);
    
    FeedFirst(t+1,1) = p.feedfirst * (1-p.muA(1)) * R(t);
    FeedFirst(t+1,2:end) = (1-p.muA(2:end-1)) .* FeedFirst(t,1:end-1);
    FeedFirstExp(t+1,2:end) = (1-p.muA(2:end-1)) .* FeedFirstExp(t,1:end-1);
    FeedFirstExp(t+1,3:4:36) = FeedFirstExp(t+1,3:4:36) + (1-mortality)  *coverage * FeedFirst(t+1,3:4:36);
    FeedFirst(t+1,3:4:36) = FeedFirst(t+1,3:4:36) - coverage * FeedFirst(t+1,3:4:36);
    
    %% calcualte the number of infectious feeder
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
%% 45% PREVALENCE - FIGURE 4B UP
%number of feeders
x=(1:t);
y1=feeders(1:t);
plot(x,y1,'-','MarkerSize',23,'LineWidth',5)
pbaspect([3 1 1]) % relative lengths of each axis.
ax = gca;
ax.FontSize = 24;% axis fontsize.
ax.LineWidth = 2;% axis LineWidth
xlabel('Time (Year)','FontSize',24) %label fontsize
ylabel('Number of Infected Mosquitoes','FontSize',24) %label fontsize
xlim([120 3285])
ylim([0 600])
set(gca,'XTick',0:365:3285)
set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9'});

%%
print -r600 -dtiff FeedersNoinsecticide_DMI_DMII_45_12X.tif

%% 45% PREVALENCE - FIGURE 4B DOWN
% infected humans
x=(1:t);
y1=IH(1:t)*100;
plot(x,y1,'-', 'MarkerSize',23,'LineWidth',4)
pbaspect([3 1 1]) % relative lengths of each axis.
ax = gca;
ax.FontSize = 24;% axis fontsize.
ax.LineWidth = 2;% axis LineWidth
xlabel('Time (Year)','FontSize',24) %label fontsize
ylabel('Infected humans (%)','FontSize',24) %label fontsize
yticks(0:10:100)
xlim([120 3285])
ylim([0 60])
set(gca,'XTick',0:365:3285)
set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9'});
%%
print -r600 -dtiff HumanNoinsecticide_DMI_DMII_45_12X.tif