
load('Model_NS3.mat')

weight_seq_all = weight_seq;
d= sum(msa_bin.*weight_seq,1)/sum(weight_seq);
load('Model_NS3_naive.mat', 'weight_seq')
load('Model_NS3_naive.mat', 'sequence_naive')
seq_bin = Binary_Seq(cell2mat(sequence_naive'),amino_single_combine_array,conserved);

d1 = sum(seq_bin .*weight_seq,1)/sum(weight_seq);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
single_MCMC = d1;
single_MSA = d;

%%
run startup.m
markersize = 3;
arrayline_min = min([single_MCMC  single_MSA]);
arrayline_max = max([single_MCMC  single_MSA]);
arrayline = arrayline_min:0.01:0.5;
FIG = figure;
color = blue;
plot(single_MSA,single_MCMC,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold;grid off;

r1 = corr(single_MSA',single_MCMC');
plot(arrayline,arrayline,'k')


xlabel('Single mutant probability (all)')
ylabel({'Single mutant probability (naive)'})
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'single_3a';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 7.84 6]);
set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
xlim([0 0.5])
ylim([0 0.5])
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%%


d= get_double_mutant(msa_bin,weight_seq_all,num_mutants_combine_array);
d1 = get_double_mutant(seq_bin,weight_seq,num_mutants_combine_array);

double_MSA_flat = d;

double_MCMC_flat = d1;
arrayline_min = min([double_MCMC_flat ;double_MSA_flat]);
arrayline_max = max([double_MCMC_flat; double_MSA_flat]);
arrayline = arrayline_min:0.01:0.5;
FIG=figure;
FIG.Name = 'double_3a';
plot(double_MSA_flat,double_MCMC_flat,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold;grid off;

r2 = corr(double_MSA_flat,double_MCMC_flat);

plot(arrayline,arrayline,'k')
xlabel('Double mutant probability (all)')
ylabel({'Double mutant probability (naive)'})
% title(['Average epsilon=' num2str(double_error)])
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 7.84 6]);
set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
xlim([0 0.5])
ylim([0 0.5])
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

%%

clear;
clc;

CON=0;
run startup.m
load('Model_NS3_naive.mat')
Single_Mutation_Oberved = sum(msa_bin.*weight_seq,1)/sum(weight_seq);
H = Single_Mutation_Oberved;
H = -log((H)./(1-H));

%% Tong2013
load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
tmp = h77;
tmp(155)='K';
seq = [h77;tmp];

tmp = h77;
tmp(36)='M';
seq = [seq;tmp];


tmp = h77;
tmp(155)='K';
tmp(36)='M';
seq = [seq;tmp];

tmp = h77;
tmp(122)='G';
seq = [seq;tmp];

tmp = h77;
tmp(122)='G';
tmp(155)='K';
seq = [seq;tmp];


tmp = h77;
tmp(383)='G';
seq = [seq;tmp];

tmp = h77;
tmp(383)='G';
tmp(115)='K';
seq = [seq;tmp];





seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        E(j) = seq_bin(j,:)*H';
    else
        E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
E = (E-mean(E))/std(E);
E1_norm = E;
I= [1;0.110000000000000;0.940000000000000;0.120000000000000;1.90000000000000;0.450000000000000;1;0.150000000000000];
% 
I1_norm=I;
I1_norm = (I1_norm-mean(I1_norm))./std(I1_norm);
[r1,p1] = corr(E1_norm,I1_norm,'type','spearman');


%% Honegger2013

load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
h77(1405-1026)='G';
seq = h77;


tmp = h77;
tmp(1403-1026)='F';
seq = [seq;tmp];


tmp = h77;
tmp(1398-1026)='R';
tmp(1409-1026)='T';
seq = [seq;tmp];




seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        E(j) = seq_bin(j,:)*H';
    else
        E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
E = (E-mean(E))/std(E);
E2_norm = E;
I= [100;17.7966101700000;20.3389830500000];
% 
I2_norm=I;
I2_norm = (I2_norm-mean(I2_norm))./std(I2_norm);
[r2,p2] = corr(E2_norm,I2_norm,'type','spearman');


%% Berger2016

load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
seq = h77;


tmp = h77;
tmp(155)='G';
seq = [seq;tmp];


tmp = h77;
tmp(155)='K';
seq = [seq;tmp];

tmp = h77;
tmp(155)='S';
seq = [seq;tmp];


tmp = h77;
tmp(168)='V';
seq = [seq;tmp];

seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        E(j) = seq_bin(j,:)*H';
    else
        E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
E = (E-mean(E))/std(E);
E3_norm = E;
I= [100;32;190;78;2.80000000000000];
% 
I3_norm=I;
I3_norm = (I3_norm-mean(I3_norm))./std(I3_norm);
[r3,p3] = corr(E3_norm,I3_norm,'type','spearman');
%% Dultz2021

load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
seq = h77;


tmp = h77;
tmp(80)='K';
seq = [seq;tmp];


tmp = h77;
tmp(80)='K';
tmp(91)='S';
tmp(174)='N';
seq = [seq;tmp];

tmp = h77;
tmp(80)='K';
tmp(174)='N';
seq = [seq;tmp];

tmp = h77;
tmp(80)='K';
tmp(91)='T';
tmp(174)='N';
seq = [seq;tmp];


seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
init_E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        init_E(j) = seq_bin(j,:)*H';
    else
        init_E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
init_I= [100.324675300000;91.5584415600000;51.6233766200000;95.4545454500000;50.6493506500000];
E=[];
I=[];
groups = [1;2;3;2;3];
for aba =1:max(groups)
        E =[E;mean(init_E(groups==aba))];
        I = [I;mean(init_I(groups==aba))];
end



% 

E = (E-mean(E))/std(E);
E4_norm = E;
I =(I-mean(I))./std(I);
I4_norm=I;
[r4,p4] = corr(E4_norm,I4_norm,'type','spearman');

%% Dultz2020_1

load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
seq =[];

tmp = h77;
tmp(41)='R';
tmp(168)='E';
seq = [seq;tmp];

tmp = h77;
tmp(168)='E';
seq = [seq;tmp];



seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
init_E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        init_E(j) = seq_bin(j,:)*H';
    else
        init_E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
init_I= [111971.831000000;62323.9436600000];
E=[];
I=[];
groups = [1;2;];
for aba =1:max(groups)
        E =[E;mean(init_E(groups==aba))];
        I = [I;mean(init_I(groups==aba))];
end



% 

E = (E-mean(E))/std(E);
E5_norm = E;
I =(I-mean(I))./std(I);
I5_norm=I;
[r5,p5] = corr(E5_norm,I5_norm,'type','spearman');

%% Dultz2020_2

load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
h77(41)='R';
seq = h77;
seq(36)='A';


tmp = h77;
tmp(36)='G';
seq = [seq;tmp];

tmp = h77;
tmp(36)='L';
seq = [seq;tmp];

tmp = h77;
tmp(36)='M';
seq = [seq;tmp];


tmp = h77;
tmp(43)='S';
seq = [seq;tmp];

tmp = h77;
tmp(54)='A';
seq = [seq;tmp];

tmp = h77;
tmp(80)='R';
seq = [seq;tmp];


tmp = h77;
tmp(138)='T';
seq = [seq;tmp];

tmp = h77;
tmp(155)='K';
seq = [seq;tmp];

tmp = h77;
tmp(155)='T';
seq = [seq;tmp];

tmp = h77;
tmp(155)='G';
seq = [seq;tmp];

tmp = h77;
tmp(155)='W';
seq = [seq;tmp];

tmp = h77;
tmp(156)='G';
seq = [seq;tmp];

tmp = h77;
tmp(156)='S';
seq = [seq;tmp];

tmp = h77;
tmp(156)='T';
seq = [seq;tmp];

tmp = h77;
tmp(156)='V';
seq = [seq;tmp];


tmp = h77;
tmp(168)='A';
seq = [seq;tmp];

tmp = h77;
tmp(168)='E';
seq = [seq;tmp];

tmp = h77;
tmp(168)='H';
seq = [seq;tmp];

tmp = h77;
tmp(168)='I';
seq = [seq;tmp];

tmp = h77;
tmp(168)='T';
seq = [seq;tmp];

tmp = h77;
tmp(168)='V';
seq = [seq;tmp];

tmp = h77;
tmp(168)='Y';
seq = [seq;tmp];

tmp = h77;
seq = [seq;tmp];

seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
init_E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        init_E(j) = seq_bin(j,:)*H';
    else
        init_E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
init_I= [50.4000000000000;2.60000000000000;76.7000000000000;27.8000000000000;3.20000000000000;14.5000000000000;24.1000000000000;1.40000000000000;23.2000000000000;4.60000000000000;2.70000000000000;2.50000000000000;41.7000000000000;62.5000000000000;0.880000000000000;1;90.5000000000000;141;80;20.1000000000000;84.7000000000000;62.5000000000000;55.4000000000000;100];
E=[];
I=[];

is_exist = logical([1;0;1;1;1;0;1;1;1;0;0;0;0;0;0;0;0;1;0;0;0;0;1;1]);
groups = [1;2;1;4;2;5;4;7;6;9;7;10;1;3;11;12;3;14;13;8;15;15;16;15];




for aba =1:max(groups)
    if ~isempty(init_E(groups==aba & is_exist))
        E =[E;mean(init_E(groups==aba & is_exist))];
        I = [I;mean(init_I(groups==aba & is_exist))];
    end
end



% 

E = (E-mean(E))/std(E);
E6_norm = E;
I =(I-mean(I))./std(I);
I6_norm=I;
[r6,p6] = corr(E6_norm,I6_norm,'type','spearman');

%% Dultz2020_2

load('data_NS3.mat', 'sequences')
h77=sequences{1,2};
seq = h77;
seq(36)='A';


tmp = h77;
tmp(36)='G';
seq = [seq;tmp];

tmp = h77;
tmp(36)='L';
seq = [seq;tmp];

tmp = h77;
tmp(36)='M';
seq = [seq;tmp];


tmp = h77;
tmp(43)='S';
seq = [seq;tmp];

tmp = h77;
tmp(54)='A';
seq = [seq;tmp];

tmp = h77;
tmp(80)='R';
seq = [seq;tmp];


tmp = h77;
tmp(138)='T';
seq = [seq;tmp];

tmp = h77;
tmp(155)='K';
seq = [seq;tmp];

tmp = h77;
tmp(155)='T';
seq = [seq;tmp];

tmp = h77;
tmp(155)='G';
seq = [seq;tmp];

tmp = h77;
tmp(155)='W';
seq = [seq;tmp];

tmp = h77;
tmp(156)='G';
seq = [seq;tmp];

tmp = h77;
tmp(156)='S';
seq = [seq;tmp];

tmp = h77;
tmp(156)='T';
seq = [seq;tmp];

tmp = h77;
tmp(156)='V';
seq = [seq;tmp];


tmp = h77;
tmp(168)='A';
seq = [seq;tmp];

tmp = h77;
tmp(168)='E';
seq = [seq;tmp];

tmp = h77;
tmp(168)='H';
seq = [seq;tmp];

tmp = h77;
tmp(168)='I';
seq = [seq;tmp];

tmp = h77;
tmp(168)='T';
seq = [seq;tmp];

tmp = h77;
tmp(168)='V';
seq = [seq;tmp];

tmp = h77;
tmp(168)='Y';
seq = [seq;tmp];

tmp = h77;
seq = [seq;tmp];

seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
init_E = zeros(size(seq,1),1);
for j =1:size(seq,1)
    if CON
        init_E(j) = seq_bin(j,:)*H';
    else
        init_E(j) = seq_bin(j,:)*J_MPF_BML*seq_bin(j,:)';

    end
end
init_I= [89.7000000000000;20.8000000000000;100.800000000000;94.7000000000000;70.2000000000000;68;51.4000000000000;0.200000000000000;47.1000000000000;71.1000000000000;45.8000000000000;28.1000000000000;45.8000000000000;50.9000000000000;2.50000000000000;0.300000000000000;34.7000000000000;97.9000000000000;55.4000000000000;43.7000000000000;97.9000000000000;47;54.2000000000000;100];
E=[];
I=[];
groups = [1;2;1;3;5;5;7;8;7;6;8;2;9;11;10;13;11;4;12;15;14;16;18;17];

is_exist = logical([1;0;1;1;1;0;1;1;1;0;0;0;0;0;0;0;0;1;0;0;0;0;1;1]);
% is_exist = ones(size(is_exist ));


for aba =1:max(groups)
    
       if  ~isempty(init_E(groups==aba & is_exist))

        E =[E;mean(init_E(groups==aba & is_exist))];
        I = [I;mean(init_I(groups==aba & is_exist))];
    end
end



% 

E = (E-mean(E))/std(E);
E7_norm = E;
I =(I-mean(I))./std(I);
I7_norm=I;
[r7,p7] = corr(E7_norm,I7_norm,'type','spearman');





%% plot figures

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

markersize = 6;
line_width = 0.25;
FIG=figure;
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
% xlabel('$${\rm E}_{\rm mutant}-{\rm E}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$${\rm f}_{\rm mutant}/{\rm f}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$$\frac{f_{mutant}}{f_{Mahoney}}$$','interpreter','latex')
hold on 
h1=plot(E1_norm,I1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',red,'LineWidth',line_width)
h2=plot(E2_norm,I2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',orange,'LineWidth',line_width)
h3=plot(E3_norm,I3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',green,'LineWidth',line_width)
h4=plot(E4_norm,I4_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',blue,'Linewidth',line_width)
h5=plot(E5_norm,I5_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
plot(E6_norm,I6_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
plot(E7_norm,I7_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
% h = legend('weigt2008','weigt2008','dunn2007','morcos2011','morcos2011');
% h.Position = [0.70 0.75 0.12 0.0488];
total_length = length(I1_norm)+length(I2_norm)+length(I3_norm)+length(I4_norm)+length(I5_norm)+length(I6_norm)+length(I7_norm)
%+length(I4_norm)+length(I5_norm)+length(I6_norm);

w1 = length(I1_norm)/total_length;
w2 = length(I2_norm)/total_length;
w3 = length(I3_norm)/total_length;
w4 = length(I4_norm)/total_length;
w5 = length(I5_norm)/total_length;
w6 = length(I6_norm)/total_length;
w7 = length(I7_norm)/total_length;
% W = [w1*ones(length(I1_norm),1);w2*ones(length(I2_norm),1);w3*ones(length(I3_norm),1);w4*ones(length(I4_norm),1);w5*ones(length(I5_norm),1);w6*ones(length(I6_norm),1)];
W = [w1*ones(length(I1_norm),1);w2*ones(length(I2_norm),1);w3*ones(length(I3_norm),1);w4*ones(length(I4_norm),1);w5*ones(length(I5_norm),1);w6*ones(length(I6_norm),1);w7*ones(length(I7_norm),1)];
rho_weighted_average = r1*w1+r2*w2+r3*w3+r4*w4+r5*w5+r6*w6+r7*w7
% +r4*w4+r5*w5+r6*w6
% text(2.1,1.5,sprintf('$$\\bar{r} = %.2f$$',rho_weighted_average),'interpreter','latex','FontSize',12)

% P = lscov([[E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm] ones(total_length,1)],[I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],W)
P = lscov([[E1_norm; E2_norm; E3_norm ; E4_norm; E5_norm; E6_norm; E7_norm ] ones(total_length,1)],[I1_norm;I2_norm;I3_norm; I4_norm;I5_norm;I6_norm;I7_norm; ],W)
% P = polyfit([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]',1);
% P = polyfit([E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm],...
%     [I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],1);
x = -2:.5:3; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1)
% h = legend('(34)','(34) M299V','(35)','(36)','(36) T292A','(37)','Location','eastoutside','Orientation','horizontal','NumColumns',1);

ylim([-2 3]);
xlim([-2 5])
xticks([-4:1:5])
% legend boxoff
if CON==1
    FIG.Name = 'conserve'; 
else
    FIG.Name = 'model'; 
end
FIG.Units = 'centimeters';
if CON==0
h = legend([h1 h2 h3 h4 h5],{'Tong2013','Honegger2013','Berger2016','Dultz2021','Dultz2020'},'Location','eastoutside','Orientation','vertical');
% h = legend('(25)','(25) M299V','(26)','(27)','(27) T292A','(28)','Location','eastoutside','Orientation','vertical');
legend boxoff
set(gcf,'Position',[10 10 10 6]);
set(gca,'Position',[.11 .18 .6 .77]);  %璋冩暣 XLABLE鍜孻LABLE涓?浼氳鍒囨�?
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
set(h,...
    'Position',[0.680104166666667 0.331023372843418 0.30646188494189 0.467344682712138],...
    'Orientation','horizontal','NumColumns',1);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
else
title("Conservation-only model")
xlabel('')
ylabel('')
ylim([-2 3]);
xlim([-2 3])
xticks([-4:1:5])
set(gcf,'Position',[10 10 4 3]);
set(gca,'Position',[.11 .18 .85 .67]);  %璋冩暣 XLABLE鍜孻LABLE涓?浼氳鍒囨�?
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

end