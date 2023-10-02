%%
run startup.m
load('Model_NS3.mat', 'sequences')
h77 = sequences{2};
load('drugs_contact.mat')
pos=[];
distance = {};
contact_site = {};
load('RAS_single.mat', 'single_ras_site')
common_ones_known = {};
common_ones_all = {};
common_ones_final = {};
final_p  =[];
for i = 1:size(seq,1)
    s = strip(seq(i,:));
    align = localalign(h77,s);
    pos = [pos;[align.Start(1),align.Stop(1)]];
    d= distance_min{1,i};   
    d = d(align.Start(2):align.Stop(2));
    distance = [distance;d];
    p = align.Start(1):align.Stop(1);
    if i==1
        align_danoprevir = align;
        p_drug_danoprevir = align.Start(2):align.Stop(2);
        p_drug_danoprevir = p_drug_danoprevir(d<=5)+981;
    end
    if i==2
        align_vaniprevir = align;
        p_drug_vaniprevir = align.Start(2):align.Stop(2);
        p_drug_vaniprevir = p_drug_vaniprevir(d<=5)+979;
    end
    if i==4
        align_grazoprevir= align;
        p_drug_grazoprevir = align.Start(2):align.Stop(2);
        p_drug_grazoprevir = p_drug_grazoprevir(d<=5)+988;
    end

    if i==5
        align_telaprevir= align;
        p_drug_telaprevir = align.Start(2):align.Stop(2);
        p_drug_telaprevir = p_drug_telaprevir(d<=5)+979;
    end
    contact_site = [contact_site;p(d<=5)];
end

strong_coupled = [36, 41, 54, 55, 71, 168,170,80,122];


%%


load('Drug_mutants.mat', 'Drugs_pos')
all_mut = Drugs_pos.Danoprevir ;
% pos = [];
% ind = [];
% mut_drug = contact_site{1,1};
% for i =1:length(all_mut)
%     p = all_mut(i);
%     if ismember(p,mut_drug)
%         ind = [ind;i];
%     end
% end
% all_mut=all_mut(ind);

all_mut = unique(all_mut);
all_p = [];


all_num = [];

% [c1,~,ia] =intersect(strong_coupled,contact_site{1,1});
interest_mut = intersect(all_mut,strong_coupled);
[c1,~,ia] =intersect(interest_mut ,contact_site{1,1});
coupled_r = p_drug_danoprevir(ia);

[c2,~,ia] =intersect(all_mut,contact_site{1,1});
all_coupled=p_drug_danoprevir(ia);
[c3,~,ia] = intersect(single_ras_site,contact_site{1,1});

% all_coupled=p_drug_danoprevir(ia);



common_ones_all = [common_ones_all;setdiff(c3,c2)];
common_ones_known = [common_ones_known; c2];
common_ones_final = [common_ones_final;c3];





p1 =pvalue_fisher(631,length(contact_site{1,1}),length(unique(all_mut)),length(intersect(all_mut,contact_site{1,1})))
n1 = length(intersect(all_mut,contact_site{1,1}));


p2 =pvalue_fisher(631,length(contact_site{1,1}),length(unique(strong_coupled)),length(intersect(strong_coupled,contact_site{1,1})))
n4 = length(intersect(strong_coupled,contact_site{1,1}));

p4 = p1;
interest_mut = intersect(all_mut,strong_coupled);
p5 = pvalue_fisher(631,length(contact_site{1,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{1,1})));
n2 = length(intersect(interest_mut,contact_site{1,1}));

p6 = pvalue_fisher(631,length(contact_site{1,1}),length(unique(single_ras_site)),length(intersect(single_ras_site,contact_site{1,1})));
n3 = length(intersect(single_ras_site,contact_site{1,1}));

p7=p2;

interest_mut = setdiff(single_ras_site,strong_coupled);
p8 = pvalue_fisher(631,length(contact_site{1,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{1,1})));

final_p = [final_p; [p4 p5 p6 p7 p8]];


rest_mut = setdiff(all_mut,strong_coupled);

p3 = pvalue_fisher(631,length(contact_site{1,1}),length(unique(rest_mut)),length(intersect(rest_mut,contact_site{1,1})))

all_p = [all_p; p1 p2 p3];
all_num = [all_num; [n1 n2 n3 n4]];
%%
load('Drug_mutants.mat', 'Drugs_pos')
all_mut = Drugs_pos.Vaniprevir ;
all_mut = unique(all_mut);

% [c1,~,ia] =intersect(strong_coupled,contact_site{2,1});
interest_mut = intersect(all_mut,strong_coupled);
[c1,~,ia] =intersect(interest_mut ,contact_site{2,1});
coupled_r2 = p_drug_vaniprevir(ia);

[c2,~,ia] =intersect(all_mut,contact_site{2,1});
all_coupled2=p_drug_vaniprevir(ia);
[c3,~,ia] = intersect(single_ras_site,contact_site{2,1});



% 
% coupled_r2 = p_drug_vaniprevir(ia);
% 
% [c2,~,ia] =intersect(all_mut,contact_site{2,1});
% all_coupled2=p_drug_vaniprevir(ia);

[c3,~,ia] = intersect(single_ras_site,contact_site{2,1});
% all_coupled2=p_drug_vaniprevir(ia);


common_ones_all = [common_ones_all;setdiff(c3,c2)];
common_ones_known = [common_ones_known; c2];
common_ones_final = [common_ones_final;c3];
p1=pvalue_fisher(631,length(contact_site{2,1}),length(unique(all_mut)),length(intersect(all_mut,contact_site{2,1})))
n1 = length(intersect(all_mut,contact_site{2,1}));



p2=pvalue_fisher(631,length(contact_site{2,1}),length(unique(strong_coupled)),length(intersect(strong_coupled,contact_site{2,1})))
n4 = length(intersect(strong_coupled,contact_site{2,1}));


p4 = p1;
interest_mut = intersect(all_mut,strong_coupled);
p5 = pvalue_fisher(631,length(contact_site{2,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{2,1})));
n2 = length(intersect(interest_mut,contact_site{2,1}));


p6 = pvalue_fisher(631,length(contact_site{2,1}),length(unique(single_ras_site)),length(intersect(single_ras_site,contact_site{2,1})));
n3 = length(intersect(single_ras_site,contact_site{2,1}));


p7=p2;
interest_mut = setdiff(single_ras_site,strong_coupled);
p8 = pvalue_fisher(631,length(contact_site{2,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{2,1})));

final_p = [final_p; [p4 p5 p6 p7 p8]];





rest_mut = setdiff(all_mut,strong_coupled);

p3 = pvalue_fisher(631,length(contact_site{2,1}),length(unique(rest_mut)),length(intersect(rest_mut,contact_site{2,1})))
all_p = [all_p; p1 p2 p3];


all_num = [all_num; [n1 n2 n3 n4]];

%%
load('Drug_mutants.mat', 'Drugs_pos')
all_mut = Drugs_pos.Telaprevir ;
all_mut = unique(all_mut);

interest_mut = intersect(all_mut,strong_coupled);
[c1,~,ia] =intersect(interest_mut ,contact_site{5,1});
coupled_r3 = p_drug_telaprevir(ia);

[c2,~,ia] =intersect(all_mut,contact_site{5,1});
all_coupled3=p_drug_telaprevir(ia);
% [c3,~,ia] = intersect(single_ras_site,contact_site{2,1});


% [c1,~,ia] =intersect(strong_coupled,contact_site{5,1});
% 
% coupled_r4 = p_drug_telaprevir(ia);
% 
% [c2,~,ia] =intersect(all_mut,contact_site{5,1});
% all_coupled4=p_drug_telaprevir(ia);
% 
[c3,~,ia] = intersect(single_ras_site,contact_site{5,1});
% all_coupled4=p_drug_telaprevir(ia);

common_ones_all = [common_ones_all;setdiff(c3,c2)];
common_ones_known = [common_ones_known; c2];
common_ones_final = [common_ones_final;c3];
p1=pvalue_fisher(631,length(contact_site{5,1}),length(unique(all_mut)),length(intersect(all_mut,contact_site{5,1})))
n1 = length(intersect(all_mut,contact_site{5,1}));



p2=pvalue_fisher(631,length(contact_site{5,1}),length(unique(strong_coupled)),length(intersect(strong_coupled,contact_site{5,1})))
n4 = length(intersect(strong_coupled,contact_site{5,1}));




p4 = p1;
interest_mut = intersect(all_mut,strong_coupled);
p5 = pvalue_fisher(631,length(contact_site{5,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{5,1})));
n2 = length(intersect(interest_mut,contact_site{5,1}));



p6 = pvalue_fisher(631,length(contact_site{5,1}),length(unique(single_ras_site)),length(intersect(single_ras_site,contact_site{5,1})));
n3 = length(intersect(single_ras_site,contact_site{5,1}));



p7=p2;
interest_mut = setdiff(single_ras_site,strong_coupled);
p8 = pvalue_fisher(631,length(contact_site{5,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{5,1})));

final_p = [final_p; [p4 p5 p6 p7 p8]];


rest_mut = setdiff(all_mut,strong_coupled);

p3 = pvalue_fisher(631,length(contact_site{5,1}),length(unique(rest_mut)),length(intersect(rest_mut,contact_site{5,1})))
all_p = [all_p; p1 p2 p3];

all_num = [all_num; [n1 n2 n3 n4]];

%%
load('Drug_mutants.mat', 'Drugs_pos')
all_mut = Drugs_pos.Grazoprevir ;
all_mut = unique(all_mut);

interest_mut = intersect(all_mut,strong_coupled);
[c1,~,ia] =intersect(interest_mut ,contact_site{4,1});
coupled_r4 = p_drug_grazoprevir(ia);

[c2,~,ia] =intersect(all_mut,contact_site{4,1});
all_coupled4=p_drug_grazoprevir(ia);

% [c1,~,ia] =intersect(strong_coupled,contact_site{4,1});
% 
% coupled_r3 = p_drug_grazoprevir(ia);
% 
% [c2,~,ia] =intersect(all_mut,contact_site{4,1});
% all_coupled3=p_drug_grazoprevir(ia);
% 
[c3,~,ia] = intersect(single_ras_site,contact_site{4,1});
% all_coupled3=p_drug_grazoprevir(ia);


common_ones_all = [common_ones_all;setdiff(c3,c2)];
common_ones_known = [common_ones_known; c2];
common_ones_final = [common_ones_final;c3];


p1=pvalue_fisher(631,length(contact_site{4,1}),length(unique(all_mut)),length(intersect(all_mut,contact_site{4,1})))
n1 = length(intersect(all_mut,contact_site{4,1}));



p2=pvalue_fisher(631,length(contact_site{4,1}),length(unique(strong_coupled)),length(intersect(strong_coupled,contact_site{4,1})))
n4 = length(intersect(strong_coupled,contact_site{4,1}));



p4 = p1;
interest_mut = intersect(all_mut,strong_coupled);
p5 = pvalue_fisher(631,length(contact_site{4,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{4,1})));
n2 = length(intersect(interest_mut,contact_site{4,1}));



p6 = pvalue_fisher(631,length(contact_site{4,1}),length(unique(single_ras_site)),length(intersect(single_ras_site,contact_site{4,1})));
n3 = length(intersect(single_ras_site,contact_site{4,1}));



p7=p2;
interest_mut = setdiff(single_ras_site,strong_coupled);
p8 = pvalue_fisher(631,length(contact_site{4,1}),length(unique(interest_mut)),length(intersect(interest_mut,contact_site{4,1})));

final_p = [final_p; [p4 p5 p6 p7 p8]];



rest_mut = setdiff(all_mut,strong_coupled);

p3 =pvalue_fisher(631,length(contact_site{4,1}),length(unique(rest_mut)),length(intersect(rest_mut,contact_site{4,1})))
all_p = [all_p; p1 p2 p3];
all_num = [all_num; [n1 n2 n3 n4]];

%%
all_p  = all_num(:,[ 1 2]);
p_values = final_p(:,[ 1 2 ]);

FIG=figure;
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
names = {'Danoprevir','Vaniprevir','Telaprevir','Grazoprevir'};


% bar(log10(1./all_p), 0.8, 'FaceColor',[0.8,0.8,0.8;0 0 0],'LineWidth',0.2,'FaceAlpha',0.4);
b = bar(all_p, 0.8);
hold on;

for i =1:size(all_p,1)
    if p_values(i,1)<0.05
        text(i-0.18,all_p(i,1)+0.2,'*','FontSize',12)
    end
    if p_values(i,2)<0.05
        text(i+0.1,all_p(i,2)+0.2,'*','FontSize',12)
    end

end



% b(1).FaceColor = [0.8,0.8,0.8];
% b(1).FaceColor = [239,118,119]/255;
% b(1).FaceColor = [0.631372549019608	0.850980392156863	0.607843137254902];
b(1).FaceColor = orange;
% b(1).FaceColor = [0.878, 0, 0];	
b(1).LineWidth = 0.2;
b(1).FaceAlpha = 0.4;

% b(2).FaceColor = [0 0 0];
% b(2).FaceColor = [0, 0, 0.886];
b(2).FaceColor  = [135,178,212]/255;
b(2).LineWidth = 0.2;
b(2).FaceAlpha = 0.4;

hold on


% plot([0 length(names)+0.5],[log10(1./0.05) log10(1./0.05)],'k--','LineWidth',0.25)
xlim([0.5 length(names)+0.5])
set(gca,'XTick',1:length(names),'XTickLabel',...
    names);

box off
hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;
set(gca,'TickDir','out')
FIG.Name = 'specific'
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Number'})
% xlabel({'Top x pairs'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';
ylim([0 15])

set(gcf,'Position',[10 10 8 6]);
set(gca,'Position',[.12 .2 .85 .77]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
legendflex([b(1),b(2)], {'DRMs','SC-DRMs'}, 'ref', gcf, ...
                       'anchor', {'ne','ne'}, ...
                       'buffer',[-10 -10], ...
                       'nrow',2, ...
                       'fontsize',8,'box','off','xscale',0.5);
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% ylim([0 4])
% set(gca,'YTick', [0:1:4])
% set(gca,'TickLength',[0.035, 0.03])

%%
all_p  = all_num(:,[ 3 4]);
p_values  = final_p(:,[ 3 4]);


FIG=figure;
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
names = {'Danoprevir','Vaniprevir','Telaprevir','Grazoprevir'};


% bar(log10(1./all_p), 0.8, 'FaceColor',[0.8,0.8,0.8;0 0 0],'LineWidth',0.2,'FaceAlpha',0.4);
b = bar(all_p, 0.8);
hold on

for i =1:size(all_p,1)
    if p_values(i,1)<0.05
        text(i-0.18,all_p(i,1)+0.2,'*','FontSize',12)
    end
    if p_values(i,2)<0.05
        text(i+0.1,all_p(i,2)+0.2,'*','FontSize',12)
    end

end


% b(1).FaceColor = [0.8,0.8,0.8];
% b(1).FaceColor = [239,118,119]/255;
b(1).FaceColor = orange;
% b(1).FaceColor = [0.631372549019608	0.850980392156863	0.607843137254902];
% b(1).FaceColor = [0.878, 0, 0];	
b(1).LineWidth = 0.2;
% b(1).FaceAlpha = 0.4;

b(2).FaceColor = blue;
% b(2).FaceColor = [0, 0, 0.886];
% b(2).FaceColor  = [135,178,212]/255;
b(2).LineWidth = 0.2;
% b(2).FaceAlpha = 0.4;

hold on


% plot([0 length(names)+0.5],[log10(1./0.05) log10(1./0.05)],'k--','LineWidth',0.25)
xlim([0.5 length(names)+0.5])
set(gca,'XTick',1:length(names),'XTickLabel',...
    names);

box off
hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;
set(gca,'TickDir','out')
FIG.Name = 'all'
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Number'})
% xlabel({'Top x pairs'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';
ylim([0 15])

set(gcf,'Position',[10 10 8 6]);
set(gca,'Position',[.12 .2 .85 .77]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
legendflex([b(1),b(2)], {'All DRMs','All SC-DRMs'}, 'ref', gcf, ...
                       'anchor', {'ne','ne'}, ...
                       'buffer',[-10 -10], ...
                       'nrow',2, ...
                       'fontsize',8,'box','off','xscale',0.5);
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% ylim([0 4])
% set(gca,'YTick', [0:1:4])
% set(gca,'TickLength',[0.035, 0.03])

