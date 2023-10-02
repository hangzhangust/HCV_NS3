run startup.m
load('ranking_couplings.mat')
load('Model_NS3.mat', 'conserved')
all_res = unique([ind_col_residue(1:300); ind_row_residue(1:300)]);
all_res = setdiff(all_res, conserved);
load('Drug_mutants.mat')
names = fieldnames(Drugs_mut);
names(ismember(names,'Ciluprevir'))=[];
names =names([6,5,7,3,9,8,2,1,4]);
names = flip(names);
names = names([9,5,8,7,4,1,3,2,6]);
p = [];
load('RAS_single.mat')

strong_coupled = [36, 41, 54, 55, 71, 168,170,80,122];
strong_coupled = setdiff(single_ras_site,strong_coupled);

for i = 1:length(names)
    Drugs_pos.(names{i}) = setdiff(Drugs_pos.(names{i}),conserved);
end


for i = 1:length(names)
p(i) = length(intersect(strong_coupled ,Drugs_pos.(names{i})));
    
    
% bar(i,p(i), 0.8, 'FaceColor',[0.2,0.2,0.2],'LineWidth',0.2,'FaceAlpha',0.6);
% hold on
end
num_coupled = p;
all_p = zeros(length(names),1);
for i = 1:length(names)
    N=515;
    n = length(all_res);
    q = length(intersect(strong_coupled ,Drugs_pos.(names{i})));
    j = length(setdiff(Drugs_pos.(names{i}), conserved));
    for k =q:min(j,n)  
        all_p(i) = all_p(i)+(nchoosek(j,k)*nchoosek(515-j,n-k))/nchoosek(N,n);
    end
end

% all_p = -log10(all_p);
% 
% all_p(1) = all_p(1)+0.1;
% 
% all_p(8) = all_p(8)-0.1;
% all_p(isinf(all_p))=0;


FIG=figure;
% plot(enrich_1a,all_p_1a,'o','MarkerEdgeColor','None','MarkerFaceColor',purple,'MarkerSize',4)
% scatter(num_coupled,all_p,'o','MarkerEdgeColor','None','MarkerFaceColor',blue,'SizeData',35);
hold on;
for i = 1:length(names)
    if all_p(i)<0.05
        text(i-0.1,num_coupled(i)+0.2,'*','FontSize',12)
    end
    
bar(i,num_coupled(i), 0.8, 'FaceColor',blue,'LineWidth',0.2,'FaceAlpha',0.6);
hold on
end
xlim([0.5 length(names)+0.5])
set(gca,'XTick',1:length(names),'XTickLabel',...
    names);
% plot(enrich_1b,all_p_1b,'o','MarkerEdgeColor','None','MarkerFaceColor',orange,'MarkerSize',5);

% plot(0:0.01:10, ones(size(0:0.01:10))*-log10(0.05), '--k','LineWidth',0.3)
% ylabel({'-log10(p-value)'})
% set(gca,'YTick', [0 1 2 3])
ylim([0 10])
% yticklabels({'0','1','2','3'})
set(gca,'TickDir','out')
set(gca,'TickLength',[0.035, 0.03])
% legend('Subtype 1a', 'Subtype 1b','Location','best')
% title('L2 reg. para. + 50%')
ylabel('Number of non-SC-DRMs')

% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'nonscdrm';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 10 8]);
set(gcf,'Position',[10 10 9.25 6.23]);
set(gca,'Position',[.15 .25 .8 .73]);  %调整 XLABLE和YLABLE不会被切掉

% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.2 .2 .78 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);