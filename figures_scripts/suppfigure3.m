%%
FIG=figure;
load('ranking_couplings.mat')

load('Drug_mutants.mat')
load('Model_NS3.mat', 'conserved')


names = fieldnames(Drugs_mut);
names(ismember(names,'ciluprevir'))=[];
names =names([6,5,7,3,9,8,2,1,4]);
names = flip(names);

for i = 1:length(names)
    Drugs_pos.(names{i}) = setdiff(Drugs_pos.(names{i}),conserved);
end
p = zeros(1,300);

for kk = 1:300
    all_p = zeros(length(names),1);
    strong_coupled = unique([ind_col_residue(1:kk); ind_row_residue(1:kk)]);
    strong_coupled = setdiff(strong_coupled, conserved);
    for i = 1:length(names)
        N=515;
        n = length(strong_coupled);
        q = length(intersect(strong_coupled ,Drugs_pos.(names{i})));
        j = length(Drugs_pos.(names{i}));
        for k =q:min(j,n)  
            all_p(i) = all_p(i)+(nchoosek(j,k)*nchoosek(515-j,n-k))/nchoosek(N,n);
        end
    end
    p(kk) = sum(all_p<=0.05);
end
hold on;



plot(p,'Color',[0.8,0.8,0.8],'LineWidth',1)
set(gca,'TickDir','out')
box off
FIG.Name = 'peak_tract'
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Number of drugs reaching' ,'the significance level'})
xlabel({'Top x pairs'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 8 6]);
set(gca,'Position',[.18 .17 .77 .8]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 9])
set(gca,'YTick', [0:3:9])
set(gca,'TickLength',[0.035, 0.03])

%%
FIG=figure;
load('ranking_couplings.mat')

load('Drug_mutants.mat')
names = fieldnames(Drugs_mut);
names(ismember(names,'ciluprevir'))=[];
names =names([6,5,7,3,9,8,2,1,4]);
names = flip(names);
for i = 1:length(names)
    Drugs_pos.(names{i}) = setdiff(Drugs_pos.(names{i}),conserved);
end


p = zeros(1,300);
for k = 1:300
    res = unique([ind_col_residue(1:k); ind_row_residue(1:k)]);
    num =0;
    for i = 1:length(names)
        if ~isempty(intersect(res ,Drugs_pos.(names{i})))
            num = num+1;
        end
    end
    p(k) = num;
end

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;




plot(p,'Color',[0.2,0.2,0.2],'LineWidth',1)
box off
set(gca,'TickDir','out')
FIG.Name = 'peak_tract'
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Number of drugs with SC-DRMs'})
xlabel({'Top x pairs'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[10 10 8 6]);
set(gca,'Position',[.18 .17 .77 .8]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')

set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
ylim([0 9])
set(gca,'YTick', [0:3:9])
set(gca,'TickLength',[0.035, 0.03])