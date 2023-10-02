r_model = [];
r_con = [];
run startup.m
for i =1:10
    filename = ['Robust_NS3/Model_' num2str(i) ".mat"];
    load(join(filename,''),"r1","r2");
    r_model = [r_model;r1];
    r_con = [r_con;r2];
end

G = [zeros(1,length(r_model)) ...
    ones(1,length(r_con))];
data = [r_model' ...
    r_con'];



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
label_xaxis_data = {'Maximum-entropy', sprintf('Conservation-only')};
text_ylabel = 'Correlation coefficient';
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

dots = r_model;
nbins = 80;
max_range = 0.8;
center = 1;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,1,50);
% f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f1=scatter(x_data,y_data ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 

dots = r_con;
nbins = 180;
max_range = 2;
center = 2;
[x_data, y_data] = dot_boxplot(dots,nbins,center,max_range,1,59);
% f2=scatter(x2,all_E(rest),'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
f2=scatter(x_data,y_data,'o','MarkerEdgeColor','w','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
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


set(gca,'YTick',-1:0.2:0)
yt = get(gca, 'YTick');
axis([xlim    -1 0])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [-0.25 -0.25], '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
P = ranksum(r_con,r_model)
% text('String','DRMs',...
%     'Position',[1.78823529411765 -68.2248677248677 0]);
text(0.85,-1.1,'model','FontSize',8);
text(1.85,-1.1,'model','FontSize',8);
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
set(gcf,'Position',[10 2  6 10]);
% set(gca,'Position',[.07 .17 .93 .8 ]);  %调整 XLABLE和YLABLE不会被切掉
set(gca,'Position',[.2 .17 .78 .8]);  %µ÷Õû XLABLEºÍYLABLE²»»á±»ÇÐµô
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.19, 0.5, 0]);
% print(['C:\Users\27909\Desktop\1'],'-dpng','-r600')
