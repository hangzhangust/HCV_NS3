load('Model_NS3.mat')
J = triu(J_MPF_BML,1);
[ ~, Ind ] = sort(J(:),1,'ascend');

[ ind_row, ind_col ] = ind2sub(size(J),Ind); % fetch indices
ia = ind_col>ind_row;
ind_col =ind_col(ia);
ind_row =ind_row(ia);

for i =1:length(ind_row)
    pos = find(num_mutants_combine_array_acc>=ind_row(i),1);
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    p = ind_row(i)-Start_site+1;
    res = amino_single_combine_array{pos,1};
    res = flip(res(2:end));
    r = res(p);
    ind_row_coverted{i} = [num2str(ind_non_conserve(pos)) r];
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
end
%%
rng default
rand_order = {};
[profile,symbols]=seqprofile(sequences);
cseq = seqconsensus(sequences);
for i =1:length(ind_non_conserve)
    p = ind_non_conserve(i);
    m = cseq(p);
    s = Symbols(find(Profile(:,p)));
    s = setdiff(s,m);
    s = s(randperm(length(s)));
    rand_order = [rand_order;s];
end
%%
ind_rand_row ={};

ind_rand_col ={};
population = 1:length(ind_non_conserve);
for i=1:1000
 
nums = randsample(population,2);
row = nums(1);
col = nums(2);
c = rand_order{row,1};
d = rand_order{col,1};
%     if isempty(c)
%         population(row)=[];
%         continue;
%     end
%     if isempty(d)
%         population(col)=[];
%         continue;
%     end
m = randperm(length(c));
n = randperm(length(d));
if ismember([num2str(ind_non_conserve(row)) c(m(1))],ind_rand_row ) && ismember([num2str(ind_non_conserve(col)) d(n(1))],ind_rand_col )
    continue;
else
    ind_rand_row = [ind_rand_row; [num2str(ind_non_conserve(row)) c(m(1))]]; 
    ind_rand_col = [ind_rand_col; [num2str(ind_non_conserve(col))  d(n(1))]];
end


%     if ~isempty(c) &&  ~isempty(d)
%         ind_rand_row = [ind_rand_row; [num2str(ind_non_conserve(row)) c(1)]];  
%         ind_rand_col = [ind_rand_col; [num2str(ind_non_conserve(col)) d(1)]];
%         c(1)=[];
%         d(1)=[];
%         rand_order{row,1}=c;
%         rand_order{col,1}=d;
%     end
end





%%
load('RAS_single.mat')
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
for i =1:length(ind_row)
    
   if  ismember(ind_col_coverted{i},single_ras) || ismember(ind_row_coverted{i},single_ras)
        count = count+1;
   end

    TPR(i) = count/i;

end




FIG=figure;
plot(TPR,'Color','k','LineWidth',1)
hold on
ind_rand_col = ind_col_coverted;
ind_rand_row = ind_row_coverted;
avg_TPR =zeros(size(ind_rand_col));
for iter =1:100
    I = randperm(length(ind_rand_col));
    ind_rand_col = ind_rand_col(I);
    I = randperm(length(ind_rand_col));
    ind_rand_row = ind_rand_row(I);
    
    TPR = zeros(size(ind_rand_col));
    count =0;
    for i =1:length(ind_rand_col)

       if  ismember(ind_rand_col{i},single_ras) || ismember(ind_rand_row{i},single_ras)
            count = count+1;
       end

        TPR(i) = count/i;

    end
    avg_TPR = avg_TPR+TPR;
end
avg_TPR = avg_TPR/100;

TPR = avg_TPR;
plot(TPR,'Color',[0.6 0.6 0.6],'LineWidth',1)

set(gca, 'XScale', 'log')
xlim([1 1e3])

xlabel('Top x pairs')
ylabel({'Precision'})
FIG.Units = 'centimeters';
FIG.Name = 'Ours';

legend('Couplings','Random')
box off
legend box off
FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
set(gcf,'Position',[10 10 10 8]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);