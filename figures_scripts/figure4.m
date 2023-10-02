%%
run startup.m
load('Model_NS3.mat')
sequences(remove_id)=[];
h77 = sequences{2};
load('Cseq_NS3.mat', 'Cseq')
Cs = h77;
Cs(ind_non_conserve) = Cseq;
Cseq = Cs;

h77(91) = 'A';

% load('RAS_single.mat', 'only_mut')
load('ranking_couplings.mat')
load('RAS_single.mat')
single_ras_converted ={};
for i =1:length(single_ras)
    if ismember(single_ras{i},only_mut)
        single_ras_converted = [single_ras_converted;single_ras{i}];
    else
     n_m = single_ras{i};   
     pos = str2num(n_m(1:end-1));
     if ismember(pos,conserved)
         continue;
     else
         m = amino_single_combine_array{ind_non_conserve==pos,1};
         m = [num2str(pos) m(end)];
         single_ras_converted =[single_ras_converted;m];
     end
    end
end
single_ras=single_ras_converted;


top_x =300;

ind_col_residue = ind_col_residue(1:top_x);
ind_row_residue = ind_row_residue(1:top_x);
ind_col_coverted = ind_col_coverted(1:top_x);
ind_row_coverted = ind_row_coverted(1:top_x);
ind_col = ind_col(1:top_x);
ind_row = ind_row(1:top_x);

mutants = {};

for i = 1:top_x
   if  ismember(ind_col_coverted{i},single_ras) 

       mutants = [mutants ;{ind_col_coverted{i} ,  ind_row_coverted{i}}];
   end
   if  ismember(ind_row_coverted{i},single_ras) 

       mutants = [mutants ;{ind_row_coverted{i} ,  ind_col_coverted{i}}];
   end
end

tmp = h77;
tmp(168) = 'E';
bin_tmp = Binary_Seq(tmp,amino_single_combine_array,conserved);
E_single = bin_tmp*J_MPF_BML*bin_tmp';
E_double_168 =[];
markers_168 = [];
mut_168 = {};
for k = 1:631
    tmp2 = tmp;
    if isempty(find(ind_non_conserve==k, 1))
%         tmp2(k) = '-';
%         bin_seq =  Binary_Seq(tmp2,amino_single_combine_array,conserved);
%         E = E_single+max(diag(J_MPF_BML));
%         E_double_168 = [E_double_168; E-E_single];
%         markers_168 = [markers_168; 0];
%         mut_168 = [mut_168; [num2str(k) '-']];
continue;
    else
        for h = setdiff(amino_single_combine_array{ind_non_conserve==k,1},tmp(k))'
            tmp2(k) = h;
            bin_seq =  Binary_Seq(tmp2,amino_single_combine_array,conserved);
            E = bin_seq*J_MPF_BML*bin_seq';
            E_double_168 = [E_double_168; E-E_single];
            mut_168 = [mut_168; [num2str(k) h]];
        if strcmp([num2str(k) h],'41X')
            markers_168 = [markers_168; 1];
        else
            markers_168 = [markers_168; 0];
        end   
        
        
        end
    end

end



tmp = h77;
tmp(80) = 'K';
bin_tmp = Binary_Seq(tmp,amino_single_combine_array,conserved);
E_single = bin_tmp*J_MPF_BML*bin_tmp';
E_double_80 =[];
markers_80 = [];
mut_80 = {};
for k = 1:631
    tmp2 = tmp;
    if isempty(find(ind_non_conserve==k, 1))
%         tmp2(k) = '-';
%         bin_seq =  Binary_Seq(tmp2,amino_single_combine_array,conserved);
%         E = E_single+max(diag(J_MPF_BML));
%         E_double_80 = [E_double_80; E-E_single];
%         markers_80 = [markers_80; 0];
%         mut_80 = [mut_80; [num2str(k) '-']];
continue;
    else
        for h = setdiff(amino_single_combine_array{ind_non_conserve==k,1},tmp(k))'
            tmp2(k) = h;
            bin_seq =  Binary_Seq(tmp2,amino_single_combine_array,conserved);
            E = bin_seq*J_MPF_BML*bin_seq';
            E_double_80 = [E_double_80; E-E_single];
            mut_80 = [mut_80; [num2str(k) h]];
        if ismember([num2str(k) h],{'91S'})
            markers_80 = [markers_80; 1];
        else
            markers_80 = [markers_80; 0];
        end   
        
        
        end
    end

end




%%


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

markersize = 6;
line_width = 0.25;
FIG=figure;

[c,b] = hist(E_double_168,50);



bar(b,c ,0.8, 'FaceColor',blue,'LineWidth',0.25,'FaceAlpha',0.6,'EdgeColor','k')
% histogram(E_double_168,'LineWidth',0.1,'BinWidth',1)
xlim([-2 15])
box off
xlabel('Energy (H77_{D168E+X}) - Energy (H77_{D168E})')
ylabel('Number')

FIG.Name = 'R1'; 
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 7.5 5]);
set(gca,'Position',[.14 .2 .83 .77]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
try 
    print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 

catch
   
end

%%


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

markersize = 6;
line_width = 0.25;
FIG=figure;

[c,b] = hist(E_double_80,50);



bar(b,c ,0.8, 'FaceColor',orange,'LineWidth',0.25,'FaceAlpha',0.6,'EdgeColor','k')
xlim([-2 15])
% histogram(E_double_80,'LineWidth',0.1,'BinWidth',1,'FaceColor',green)
box off
xlabel('Energy (H77_{Q80K+X}) - Energy (H77_{Q80K})')
ylabel('Number')

FIG.Name = 'R2'; 
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 7.5 5]);
set(gca,'Position',[.14 .2 .83 .77]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
try 
    print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 

catch
   
end
