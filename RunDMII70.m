%% RunDMII70.m
%% Required functions:
% DM.m
% DMParameters.m
% mortality_ot.m
% brewermap.m
%%
clc;clear;
%% 
p = DMParameters();
p.cov_vec = 0:.025:1;
IR_vec = 0:.025:1;
%%
mortality_out = zeros(1, 41);
DMII = NaN(length(IR_vec),length(p.cov_vec));
DMI = NaN(length(IR_vec),length(p.cov_vec));

p.bites = 152/p.K; 
disp(p.bites)
%%
for iterate = 1:length(IR_vec)
    
    p.insect = IR_vec(iterate);
    disp(p.insect)
    
    p.hs = 0;
    [IHout,cov_vec] = DM(p);
    DMI(iterate,:) = (IHout(1)-IHout)/IHout(1);
    
    IR=IR_vec(iterate);
    mortality_out(iterate)=mortality_ot(IR);
    
    p.insect = mortality_out(iterate); % p.insect is the mortality of insects on contact with insecticides
    [IHout,cov_vec] = DM(p);
    DMII(iterate,:) = (IHout(1)-IHout)/IHout(1);
end

%% 70% PREVALENCE - FIGURE S9 LEFT

figure()
colormap(brewermap([],'YlGnBu'))
imagesc(cov_vec,1-IR_vec,flipud(DMI))
xlabel('Coverage (%)')
ylabel('Insecticide Resistance (%)')
set(gca,'fontsize',18,'fontweight','normal','xtick',0:.2:1,'xticklabel',0:20:100,'ytick',0:.2:1,'yticklabel',fliplr(0:20:100))
axis square
title('Deltamethrin I','fontweight','normal')
%%
print -r600 -dtiff 12172019_1_70_12X_YlGnBu_1.tif

%% 70% PREVALENCE - FIGURE S9 RIGHT

figure()
colormap(brewermap([],'YlGnBu'))
imagesc(cov_vec,1-IR_vec,flipud(DMII))
xlabel('Coverage (%)')
ylabel('Insecticide Resistance (%)')
set(gca,'fontsize',18,'fontweight','normal','xtick',0:.2:1,'xticklabel',0:20:100,'ytick',0:.2:1,'yticklabel',fliplr(0:20:100))
c= colorbar;
c.Label.String = 'Effectiveness (%)';
c.Ticks = 0:.1:1;
c.TickLabels = 0:10:100;
c.FontSize  = 18;
title('Deltamethrin II','fontweight','normal')
axis square
%%
print -r600 -dtiff 12172019_2_70_12X_YlGnBu_1.tif
