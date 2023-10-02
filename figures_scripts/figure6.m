%%
run startup.m
load('Drug_mutants.mat')
load('mean_escape_time_NS3.mat')
% load('delatE_avg_min.mat')
all_E_avg = mean_escape_time_NS3;
load('RAS_single.mat', 'single_ras_site')
load('Model_NS3.mat', 'conserved')
single_ras_site = setdiff(single_ras_site, conserved);
delta_E_RAS = all_E_avg(single_ras_site);
[delta_E_RAS,I] = sort(delta_E_RAS,'ascend');
names= {};
for i = 1:length(single_ras_site)
    names = [names;num2str(single_ras_site(i))];
end


names  =names(I);
single_ras_site = single_ras_site(I);
FIG=figure;


hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;
% load('color_gradient.mat')
names_compen = [36,41,54,55,71,80,122,168,170];
% for i =1:length(names)
% c = min(floor(delta_E_RAS(i)/200*size(gradient,1)),100);
% bar(i,delta_E_RAS(i), 0.5, 'FaceColor',gradient(c,:),'LineWidth',0.2,'FaceAlpha',0.8);
% 
% hold on
% 
% end
for i =1:length(names)
% c = min(floor(delta_E_RAS(i)/200*size(gradient,1)),100);
% bar(i,delta_E_RAS(i), 0.5, 'FaceColor',gradient(c,:),'LineWidth',0.2,'FaceAlpha',0.8);
if ismember(single_ras_site(i),names_compen )
     b1 = bar(i,delta_E_RAS(i), 0.5, 'FaceColor',blue,'LineWidth',0.2,'FaceAlpha',0.6);
else
    b2 = bar(i,delta_E_RAS(i), 0.5, 'FaceColor',red,'LineWidth',0.2,'FaceAlpha',0.6);
end
hold on

end



% plot([0 length(names)+1],[100 100],'k--')
xlim([0.5 length(names)+0.5])
set(gca,'XTick',1:length(names),'XTickLabel',names);
set(gca,'TickDir','out')
FIG.Name = 'deltaE';
% set(gca,'TickLength',[0.02, 0.01])
% ylabel({'Escape time'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';

% xtickangle(45)
set(gcf,'Position',[10 2  13 5]);
set(gca,'Position',[.07 .17 .93 .8 ]);  %调整 XLABLE和YLABLE不会被切掉

legendflex([b1,b2], {'SC-DRMs','Remaining DRMs'}, 'ref', gcf, ...
                       'anchor', {'nw','nw'}, ...
                       'buffer',[50 -10], ...
                       'nrow',2, ...
                       'fontsize',8,'box','off','xscale',0.5);

% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
% ytickangle(45)
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.05, 0.5,0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% ylim([0 17])
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
%%

run startup.m
rng default

load('Model_NS3.mat')
load('RAS_single.mat')
load('mean_escape_time_NS3.mat')
deltaE = mean_escape_time_NS3;
L=length(deltaE);
all_E=mean_escape_time_NS3;


load('Drug_mutants.mat')
run startup.m
load('RAS_single.mat', 'single_ras_site')
load('Model_NS3.mat', 'conserved')
single_ras_site = setdiff(single_ras_site, conserved);
polymorphisms_associated_with_neutralization_resistance =[36,41,54,55,71,80,122,168,170];

rest = setdiff(single_ras_site,polymorphisms_associated_with_neutralization_resistance);
G = [zeros(1,length(polymorphisms_associated_with_neutralization_resistance)) ...
    ones(1,length(rest))];
data = [all_E(polymorphisms_associated_with_neutralization_resistance) ...
    all_E(rest)];

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
L=363;
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'SC-DRMs', sprintf('Remaining')};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 500;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')
size_marker = 20;
x1=0.8+0.4*(rand(length(polymorphisms_associated_with_neutralization_resistance),1));
x2=1.8+0.4*(rand(length(rest),1));
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 


dots = all_E(polymorphisms_associated_with_neutralization_resistance);
nbins = 80;
max_range = 0.8;
center = 1;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,1,50);
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f1=scatter(x_data,y_data ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 

dots = all_E(rest);
nbins = 180;
max_range = 2;
center = 2;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,1,59);
% f2=scatter(x2,all_E(rest),'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
f2=scatter(x_data,y_data,'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
hold on
figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);
hold on;


% red = color_scheme_npg(1,:);
% blue = color_scheme_npg(2,:);
% V=violinplot(data, [repmat("Escape",1,length(polymorphisms_associated_with_neutralization_resistance)) repmat("Remaining",1,length(rest))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel('Escape time');
% V(1, 1).ViolinColor = blue;
% V(1, 2).ViolinColor = red;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;
% V(1, 1).BoxColor  = 'None';
% V(1, 1).BoxColor  = 'None';
set(gca,'YTick',0:75:300)
yt = get(gca, 'YTick');
axis([xlim    0  400])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*330, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
P = ranksum(all_E(rest),all_E(polymorphisms_associated_with_neutralization_resistance),'tail','right')
P = ranksum(all_E(rest),all_E(polymorphisms_associated_with_neutralization_resistance))
text('String','DRMs',...
    'Position',[1.78823529411765 -68.2248677248677 0]);
% text(0.824561403508772,-67.46842105263158,'RAS','FontSize',8);
% ind = floor(log10(P));
% P = roundn(P/10^ind,-1);
% t = ['$$ P = ',num2str(P),' \times 10^{',num2str(ind),'} $$'];
% plot(xt([1 2]), [1 1]*max(yt), '-k','LineWidth',0.5)
% text(1,500,t,'interpreter','latex','FontSize',8)
% text(1,max(yt)*1.15,'$$ P = 9.3 \times 10^{-20} $$','interpreter','latex','FontSize',8)
% text(2,-100,'DRMs','FontSize',8,'FontName', 'Arial')

set(gca,'TickLength',[0.02, 0.01])
FIG.Name = 'Escape_mutation';
set(gca,'TickDir','out')
FIG.Units = 'centimeters';
FIG.Name = 'Escape_compare';
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gcf,'Position',[10 2  5 5]);
% set(gca,'Position',[.07 .17 .93 .8 ]);  %调整 XLABLE和YLABLE不会被切掉
set(gca,'Position',[.25 .17 .72 .8]);  %µ÷Õû XLABLEºÍYLABLE²»»á±»ÇÐµô
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.235, 0.5, 0]);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end

% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
% set(gca,'FontName','Arial','FontSize',8)


% % export_fig C:\Users\27909\Desktop\Escape_mutation.png -native