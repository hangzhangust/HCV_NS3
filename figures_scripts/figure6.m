load('Drug_muants.mat')
load('mean_escape_time_NS3.mat')
% load('delatE_avg_pre_filter.mat')
% Drugs_SVR = Drugs_SVR_multi;
% f = fieldnames(Drugs_SVR_multi);
%  for i = 1:length(f)
%     Drugs_SVR.(f{i}) = Drugs_SVR_multi.(f{i});
%  end
deltaE = mean_escape_time_NS3_single;
run startup.m
load('Model_NS3.mat')
names = fieldnames(Drugs_mut);
names(ismember(names,'ciluprevir'))=[];
names =names([6,5,7,3,9,8,2,1,4]);
names = flip(names);
data_drugs={};
all_mut = [];
for n =names'
    mut = Drugs_pos.(n{1});
    
     all_mut = [all_mut;mut];
     data_drugs = [data_drugs mut];
end


e = [];
tc = [36,41,54,55,71,80,122,168,170];
% load('RAS_single.mat', 'single_ras_site')
% tc = single_ras_site;
for kk = 1:length(data_drugs)
    p = length(intersect(data_drugs{kk},tc));
    e = [e p];
end
% for kk = 1:length(data_drugs)
%     p = intersect(data_drugs{kk},tc);
%     e = [e max(mean_escape_time_NS3(p))];
% end



svr=[];
std_svr = [];
[~,ia,ib] = intersect(fieldnames(Drugs_SVR),names);
names = names(ib);
e = e(ib);
for n =names'
         d = Drugs_SVR.(n{1});
         s = d(:,1);
         weight = d(:,2);
         std_svr = [std_svr; std(s,weight)]; 
         svr =[svr;sum(s.*weight)/sum(weight)]; 
%          svr =[svr;mean(s)]; 
end
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
FIG=figure;

P = lscov([e' ones(length(svr),1)],svr)
% P = polyfit([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]',1);
% P = polyfit([E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm],...
%     [I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],1);
x = 2:1:6; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1)
hold on
markersize=5;
ylim([60 80])
yticks([50:10:90])
xlim([2 6])
xticks([2:6])
color=[0.6 0.6 0.6];
[r,pval]=corr(e',svr,'type','spearman')
mean_escape_time_E1E2 =e;
mean_escape_time=svr';
% arrayline_min = min([mean_escape_time_E1E2 mean_escape_time]);
% arrayline_max = max([mean_escape_time_E1E2 mean_escape_time]);
% arrayline = arrayline_min:0.01:arrayline_max;
plot(mean_escape_time_E1E2,mean_escape_time,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold;grid off;
% plot(arrayline,arrayline,'k')

hold on;

% 
% err= errorbar(mean_escape_time_E1E2,mean_escape_time, std_svr, 'LineStyle','none');
% 
% err.Color = [0 0 0];
% err.LineWidth = 0.5;



FIG.Name = 'Escape_corr';
set(gca,'TickDir','out')
FIG.Units = 'centimeters';
FIG.Name = 'crr';
set(gcf,'Position',[6.53 6.53 7.84 6]);
box off
% text(6.5,61,sprintf('r = %.2f',r),'FontSize',12)
xlabel({'Number of SC-DRMs'})
ylabel({'Sustained virological response rates (%)'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.15, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');

%%
load('mean_escape_time_NS3.mat')
load('Drug_muants.mat')
% load('deltaE_avg.mat')
Drugs_SVR = Drugs_SVR_multi;
% f = fieldnames(Drugs_SVR_multi);
%  for i = 1:length(f)
%     Drugs_SVR.(f{i}) = Drugs_SVR_multi.(f{i});
%  end
% deltaE = deltaE_avg;
deltaE = mean_escape_time_NS3;
run startup.m
load('Model_NS3.mat')
names = fieldnames(Drugs_mut);
names(ismember(names,'ciluprevir'))=[];
names =names([6,5,7,3,9,8,2,1,4]);
names = flip(names);
% names(ismember(names,'simeprevir'))=[];
% names(ismember(names,'danoprevir'))=[];
% names(ismember(names,'grazoprevir'))=[];
% names(ismember(names,'glecaprevir'))=[];


data_drugs={};
all_mut = {};
for n =names'
    mut = Drugs_pos.(n{1});

     
     
     data_drugs = [data_drugs; mut];
end


e = [];
tc = [36,41,54,55,71,80,122,168,170];
% load('RAS_single.mat', 'single_ras_site')
% tc = single_ras_site;
for kk = 1:length(data_drugs)
    p = length(intersect(data_drugs{kk},tc));
    e = [e p];
end


svr=[];
[~,ia,ib] = intersect(fieldnames(Drugs_SVR),names);
names = names(ib);
e = e(ib);
std_svr = [];
for n =names'
         d = Drugs_SVR.(n{1});
         s = d(:,1);
         weight = d(:,2);
         svr =[svr;sum(s.*weight)/sum(weight)]; 
         std_svr = [std_svr; std(s,weight)]; 
%          svr =[svr;mean(s)]; 
end
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
FIG=figure;

P = lscov([e' ones(length(svr),1)],svr)
% P = polyfit([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]',1);
% P = polyfit([E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm],...
%     [I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],1);
x = 2:1:8; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1)
hold on
markersize=5;
ylim([90 100])
yticks([85:5:100])
ylim([85 100])
% yticks([90:5:100])
% xlim([50 200])
color=[0.6 0.6 0.6];
[r,pval]=corr(e',svr,'type','spearman')
mean_escape_time_E1E2 =e;
mean_escape_time=svr';
% arrayline_min = min([mean_escape_time_E1E2 mean_escape_time]);
% arrayline_max = max([mean_escape_time_E1E2 mean_escape_time]);
% arrayline = arrayline_min:0.01:arrayline_max;
plot(mean_escape_time_E1E2,mean_escape_time,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold;grid off;
% plot(arrayline,arrayline,'k')
hold on
% err= errorbar(mean_escape_time_E1E2,mean_escape_time, std_svr, 'LineStyle','none');
% 
% err.Color = [0 0 0];
% err.LineWidth = 0.5;

FIG.Name = 'Escape_corr';
set(gca,'TickDir','out')
FIG.Units = 'centimeters';
FIG.Name = 'multi';
% xlim([4 9])
set(gcf,'Position',[6.53 6.53 5 10]);
box off
text(6.5,61,sprintf('r = %.2f',r),'FontSize',12)
set(gcf,'Position',[6.53 6.53 7.84 6]);
box off
% text(6.5,61,sprintf('r = %.2f',r),'FontSize',12)
xlabel({'Number of SC-DRMs'})
ylabel({'Sustained virological response rates (%)'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.15, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');