load('Model_NS3.mat')
J = triu(J_MPF_BML,1);
[ ~, Ind ] = sort(J(:),1,'ascend');

[ ind_row, ind_col ] = ind2sub(size(J),Ind); % fetch indices
ia = ind_col>ind_row;
ind_col =ind_col(ia);
ind_row =ind_row(ia);
all_res=[];
ind_col_residue=zeros(length(ind_row),1);
ind_row_residue=zeros(length(ind_row),1);
for i =1:length(ind_row)
    pos = find(num_mutants_combine_array_acc>=ind_row(i),1);
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    p = ind_row(i)-Start_site+1;
    res = amino_single_combine_array{pos,1};
    res = flip(res(2:end));
    r = res(p);
    ind_row_coverted{i} = [num2str(ind_non_conserve(pos)) r];
    ind_row_residue(i) = ind_non_conserve(pos);
end
for i =1:length(ind_col)
    pos = find(num_mutants_combine_array_acc>=ind_col(i),1);
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    p = ind_col(i)-Start_site+1;
    res = amino_single_combine_array{pos,1};
    res = flip(res(2:end));
    r = res(p);
    ind_col_coverted{i} = [num2str(ind_non_conserve(pos)) r];
    ind_col_residue(i) = ind_non_conserve(pos);
end

%%
run startup.m
load('RAS_single.mat')
% load('deltaE_avg.mat', 'only_mut')
TPR = zeros(size(ind_row));
single_ras_conerted ={};
for i =1:length(single_ras)
    if ismember(single_ras{i},only_mut)
        single_ras_conerted = [single_ras_conerted;single_ras{i}];
    else
     n_m = single_ras{i};
     pos = str2num(n_m(1:end-1));
     if ismember(pos,conserved)
         continue;
     else
         m = amino_single_combine_array{ind_non_conserve==pos,1};
         m = [num2str(pos) m(end)];
         single_ras_conerted =[single_ras_conerted;m];
     end
    end
end
single_ras=single_ras_conerted;

count =0;
all_p=[];
all_n = [];

for nums=[10 100 300]
markersize = 6;
line_width = 2;

names={};
for i =1:631
    names = [names;num2str(i)];
end
angles =[];
for i =1:length(names)
    angles = [angles;(i-1)/length(names)*360];
end
pos_ind=[];
neg_ind =[];
name_ind=[];

all_res = [];
for i =1:nums
       resi = only_mut(ind_row(i));
       resi = resi{1};
   resi = resi(1:end-1);
   all_res = [all_res;str2num(resi)];
   resj = only_mut(ind_col(i));
   resj = resj{1};
   resj = resj(1:end-1);
   all_res = [all_res;str2num(resj)];
   if  ismember(ind_col_coverted{i},single_ras) || ismember(ind_row_coverted{i},single_ras)

       pos_ind = [pos_ind ;[find(strcmp(resi,names))    find(strcmp(resj,names))]];
   else   
       neg_ind = [neg_ind ;[find(strcmp(resi,names))    find(strcmp(resj,names))]];
   end
   if ismember(ind_col_coverted{i},single_ras)
       name_ind = [name_ind ; find(strcmp(resj,names))];
   end
   if ismember(ind_row_coverted{i},single_ras)
       name_ind = [name_ind ; find(strcmp(resi,names))];
   end

end
r=1.1;
x= r*sin([0:1:360]/360*2*pi);
y = r*cos([0:1:360]/360*2*pi);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)



    
FIG=figure;
plot(x,y,'Color',[0 0 0]);
hold on;
for i = 1:size(neg_ind,1)
    alpha = angles(neg_ind(i,1));
    beta = angles(neg_ind(i,2));
    color = [0.8 0.8 0.8];
    if abs(alpha-beta)~=180
        plot_arc(alpha,beta,color,line_width)
    else
        u  = [cos(alpha/180*pi),cos(beta/180*pi)];
        v  = [sin(alpha/180*pi);sin(beta/180*pi)];
        plot(u,v,'Color',color,'LineWidth', line_width)
    end
    
end
for i = 1:size(pos_ind,1)
    alpha = angles(pos_ind(i,1));
    beta = angles(pos_ind(i,2));
    if abs(alpha-beta)~=180
    plot_arc(alpha,beta,orange,line_width)
    else
        u  = [cos(alpha/180*pi);sin(alpha/180*pi)];
        v  = [cos(beta/180*pi);sin(beta/180*pi)];
        plot(u,v,'Color',color,'LineWidth', line_width)
    end    
    
end
name_ind = unique(name_ind);
for i=1:length(name_ind)
    ang = angles(name_ind(i));
    v  = 1.05*r*[cos(ang/180*pi);sin(ang/180*pi)];
    align ='left';
    if strcmp(names{name_ind(i)},'168') && any(contains(names(name_ind),'170'))
        txt=text(v(1)+0.037,v(2)+0.01,names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
        continue
    end
    if strcmp(names{name_ind(i)},'170') && any(contains(names(name_ind),'168'))
        txt=text(v(1)-0.037,v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
        continue
    end
    if strcmp(names{name_ind(i)},'55') && any(contains(names(name_ind),'54'))
        txt=text(v(1)-0.03,v(2)+0.03,names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
        continue
    end
    if strcmp(names{name_ind(i)},'54') && any(contains(names(name_ind),'55'))
        txt=text(v(1)+0.03,v(2)-0.03,names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
        continue
    end
    if strcmp(names{name_ind(i)},'36') && any(contains(names(name_ind),'41'))
        txt=text(v(1)+0.015,v(2)-0.0275,names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
        continue
    end
%     if  angles(i)>90 && angles(i) <270
%         angles(i) = angles(i)+180;
%         align ='right';
%     end
%     angles(i) = mod(angles(i),180);
    txt=text(v(1),v(2),names{name_ind(i)},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',ang);
end

pos_res =[];
for i =1:length(name_ind)
    pos_res = [pos_res;str2num(names{name_ind(i)})];
end
pos_res = unique(pos_res);


set(gca,'XColor', 'none','YColor','none')
xlim([-1.4,1.4])
ylim([-1.4,1.4])
FIG.Name = join(["Top_",num2str(nums)],'');
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 7 7]);
set(gca,'Position',[.01 .01 .98 .98]);  %è°æ´ XLABLEåYLABLEä¸?ä¼è¢«åæ?
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])

% try 
%     print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 
% 
% catch
%     try
%         print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
%     catch 
%         print(['/Users/mac/Documents/' FIG.Name],'-dpng','-r600');
%     end
%         
% end



% for i=1:length(angles)
%     v  = 1.05*r*[cos(angles(i)/180*pi);sin(angles(i)/180*pi)];
%     align ='left';
% %     if  angles(i)>90 && angles(i) <270
% %         angles(i) = angles(i)+180;
% %         align ='right';
% %     end
% %     angles(i) = mod(angles(i),180);
%     txt=text(v(1),v(2),names{i},'HorizontalAlignment',align,'VerticalAlignment', 'middle','rotation',angles(i));
% end


load('Drug_mutants.mat');
names = fieldnames(Drugs_pos);
names(ismember(names,'ciluprevir'))=[];
names =names([6,5,7,3,9,8,2,1,4]);

names = flip(names);
info={};
p_values = [];
num_mut = [];
for n =names'
    mut = Drugs_pos.(n{1});
         p=pvalue_fisher(631,length(unique(all_res)),length(unique(mut)),length(intersect(mut,pos_res)));
         info = [info; n p length(intersect(mut,pos_res))];
         p_values = [p_values; p ];
         num_mut = [num_mut; length(intersect(mut,pos_res))];

end
all_p = [all_p p_values];
all_n = [all_n num_mut];
end


%% Sectors;
clear
run startup.m
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultTextFontSize',6)

N = 631;
num_PC='7';
load('RAS_single.mat', 'single_ras_site')
name = ['RocaSecs_numPCs' num_PC '.csv'];

T=csvread(name);
inter_ras={};
for i =1:str2num(num_PC)
    tmp = T(i,:);
    tmp = tmp(tmp~=0);
    sec{i} = tmp;
    inter_ras{i} = intersect(tmp,single_ras_site);
end




if size(sec,1)==1
    sec=sec';
end
% sec{2} = [sec{1,1} sec{2,1}];
% sec{4} = [sec{4,1} sec{5,1} sec{6,1} sec{7,1}];
% sec = sec(2:4);


region = [];
region{1,1} = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22];
region{2,1} = 10:24;
region{3,1} = [1:10 139];
region{4,1} = [460:467];
region{5,1} = [435;436;437;438;439;440;441;442;443;444;445;446;447;448;449;450;451;452;453;477;478;479;480;481];
names_compen = [36,41,54,55,71,80,122,168,170];

region{6,1} = single_ras_site;
region{7,1} = names_compen;
p_value = zeros(length(sec),length(region));

site_in_each = {};

for j =1:length(sec)
for i = 1:length(region)
    	
    x = [length(intersect(sec{j,1},region{i,1}));length(intersect(setdiff(1:N,sec{j,1}),region{i,1}))];
    y = [length(intersect(sec{j,1},setdiff(1:N,region{i,1})));length(intersect(setdiff(1:N,sec{j,1}),setdiff(1:N,region{i,1})))];
    tmp = table(x,y);
    [h,p,stats] = fishertest(tmp);
     p_value(j,i)=pvalue_fisher(631,length(sec{j,1}),length(region{i,1}),length(intersect(sec{j,1},region{i,1})));
%     p_value(j,i) = p;

   if i==7 
       if length(intersect(region{i,1},sec{j,1}))>0
site_in_each = [site_in_each; intersect(region{i,1},sec{j,1})];
       else
           site_in_each = [site_in_each; [1]];
       end
   end
end
end

FIG=figure;



% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end

% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
% set(gca,'FontName','Arial','FontSize',8)


% % export_fig C:\Users\27909\Desktop\Escape_mutation.png -native
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');

a = brewermap(100,'Blues').';
b = fliplr(a(:,1:length(a)));
c = b.';
color_scheme_PuBuGnFlipped = c;
label_yaxis_data={};
for i=1:length(sec)
    label_yaxis_data = [label_yaxis_data;num2str(i)];
end
label_xaxis_data={'NS3-NS4A-Pro-Act','NS3-NS4A-Mem-Asso','NS5A-Hyper-Phos','NS3-Motif-Enz-Heli','NS3-Intra-Dimer-Int','DRMs','SC-DRMs'};
% label_xaxis_data={'NS3-NS4A-Pro-Act','NS3-NS4A-Mem-Asso','NS5A-Hyper-Phos','NS3-Motif-Enz-Heli','NS3-Intra-Dimer-Int'};

h=heatmap(label_xaxis_data,label_yaxis_data,p_value,'GridVisible','on','ColorLimits',[0 0.05],'ColorbarVisible','off','FontName','Arial','FontSize',7,'Colormap',color_scheme_PuBuGnFlipped);
FIG.Name = 'heatmap';
FIG.Units = 'centimeters';
FIG.Name = 'crr';
set(gcf,'Position',[6.53 6.53 10 8]);
% text(6.5,61,sprintf('r = %.2f',r),'FontSize',12
ylabel({'Sector'})
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'Position',[.13 .19 .86 .8]);  %调整 XLABLE和YLABLE不会被切掉
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.17, 0.5, 0]);

% set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);