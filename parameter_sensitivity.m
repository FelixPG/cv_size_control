%% This script computes the numerical Hessian of a particular measure of your choice (e.g. CV or mean)
% and diagonalizes the matrix to find the eigenvalues and eigenvectors,
% corresponding to parameter sensitivity directions. Since only 1 point is
% differentiated in the matrix, the first eigenvector corresponds to the
% gradient of the measure in parameter space. 

close all
clear all

parameters = {'start_trans_slope';'start_trans_thresh';'PSG1_slope';'PSG1_thresh';'div_trans_d_slope_1';...
         'div_trans_d_slope_2';'div_trans_d_thresh';'growth_rate';'PSG1_mother_slope';'PSG1_mother_thresh';...
         'div_trans_m_slope_1';'div_trans_m_slope_2';'div_trans_m_thresh';'m_growth_rate_1';'m_growth_rate_2';'d_bud_mass_correction';
         'sg2m_timer_slope';'sg2m_timer_thresh';'m_g2_timer_slope';'m_g2_timer_thresh';'m_bud_mass_correction'};

parameter_values = [1.85056126550248e-06; 0.0337256466753454; -2.596953686853002e-04; 34.316569344194875; 7.628399268975367e-06; ...
         -2.983508662043306e-06; -0.011854207243603; 0.0049254433; -5.429621092521212e-04; 47.735552938244965; ...
            6.610122886272357e-06; -2.538505485426652e-06; -0.006267310306677; 0.001953883911940; 98.154228802411230; 1.369421416415018e+03;
            -2.681657833310124e-04; 1.094494905033565e+02; -2.979402851076718e-04; 1.049573409490653e+02; 2.417144252429131e+03];

row_labels = {'START slope(M)','START offset','Post-START G1 (daughter) slope (M)','Post-START G1 (daughter) offset','Division (daughter) slope (M_{bud})',...
    'Division (daughter) slope (M_{START})','Division (daughter) offset','Constant growth rate','Post-START G1 (mother) slope (M)','Post-START G1 (mother) offset',...
    'Division (mother) slope (M_{bud})','Division (mother) slope (M_{START})','Division (mother) offset','Growth rate (mother) 1','Growth rate (mother) 2',...
    'Bud mass correction (daughter)','S/G2/M timer (daughter) slope (M_{START})',...
    'S/G2/M timer (daughter) offset','S/G2/M timer (mother) slope (M_{START})','S/G2/M timer (mother) offset','Bud mass correction (mother)'};

mean_size_birth_wt = 2.737635054332483e+04;
cv_size_birth_wt = 0.301088407778900;
std_size_birth_wt = 8.242701795886689e+03;

mean_size_start_wt = 3.392439354741067e+04;
cv_size_start_wt = 0.243939728008413;
std_size_start_wt = 8.275507334805718e+03;

mean_size_cyto_wt = 6.212302878327006e+04;
cv_size_cyto_wt = 0.219461070776178;
std_size_cyto_wt = 1.363358641663578e+04;

mean_size_full_wt = 4.770394548246391e+04;
cv_size_full_wt = 0.415597997316403;
std_size_full_wt = 1.982566420660287e+04;


measure = {'mean';'std';'cv'};

%%%%%%%%%%%%%%%%%%%%
% Indicate here which distribution you would like to study. It will look at
% the mean, the variance and the CV at the same time. 

which_dist = 'full'; % Choices include 'full','birth','start','cyto'

%%%%%%%%%%%%%%%%%%%%

bins = -0.5:0.025:0.5;%0.5e+04:0.25e+04:12.5e+04;
num_p = length(parameter_values);
m = zeros(num_p,num_p,3);
jac = zeros(num_p,3);

% data_path corresponds to the point around which you would like to
% numerically differentiate. In principle, this should be the where the
% results of the baseline WT simulation should be.
data_path = '.../';
% Change save path as needed
save_path = '.../';
%% 
format long
for q = 1:num_p
    path_q = strcat(data_path,char(parameters(q)));
    load(strcat(path_q,'/ph/derivatives_workspace.mat'),'cell_v','cell_bud','size_at_birth','volume_at_START','cell_cycle_growth');

    if strcmp(which_dist,'full')
        size_q_p = cell_v(:,length(cell_v(1,:))) + cell_bud(:,length(cell_v(1,:)));
        mean_wt = mean_size_full_wt;
        std_wt = std_size_full_wt;
        cv_wt = cv_size_full_wt;
    elseif strcmp(which_dist,'birth')
        size_q_p = size_at_birth;
        mean_wt = mean_size_birth_wt;
        std_wt = std_size_birth_wt;
        cv_wt = cv_size_birth_wt;
    elseif strcmp(which_dist,'start')
        size_q_p = volume_at_START;
        mean_wt = mean_size_start_wt;
        std_wt = std_size_start_wt;
        cv_wt = cv_size_start_wt;
    elseif strcmp(which_dist,'cyto')
        size_q_p = size_at_birth + cell_cycle_growth;
        mean_wt = mean_size_cyto_wt;
        std_wt = std_size_cyto_wt;
        cv_wt = cv_size_cyto_wt;
    end
    
    mean_q_p = mean(size_q_p);%;0.5*(mean(size_q_p) - mean_wt)^2
    std_q_p = std(size_q_p);%;0.5*(std(size_q_p) - std_wt)^2
    cv_q_p = std_q_p/mean_q_p;%(std_q_p/mean_q_p - cv_wt)^2;

    load(strcat(path_q,'/mh/derivatives_workspace.mat'),'cell_v','cell_bud','size_at_birth','volume_at_START','cell_cycle_growth');
    
    if strcmp(which_dist,'full')
        size_q_m = cell_v(:,length(cell_v(1,:))) + cell_bud(:,length(cell_v(1,:)));
    elseif strcmp(which_dist,'birth')
        size_q_m = size_at_birth;
    elseif strcmp(which_dist,'start')
        size_q_m = volume_at_START;
    elseif strcmp(which_dist,'cyto')
        size_q_m = size_at_birth + cell_cycle_growth;
    end

    mean_q_m = mean(size_q_m);%(mean(size_q_m) - mean_wt)^2;
    std_q_m = std(size_q_m);%(std(size_q_m) - std_wt)^2;
    cv_q_m = std_q_m/mean_q_m;%(std_q_m/mean_q_m - cv_wt)^2;

    jac(q,1) = (mean_q_p - mean_q_m)/(2*0.05);%derivative wrt log(param)
    jac(q,2) = (std_q_p - std_q_m)/(2*0.05);
    jac(q,3) = (cv_q_p - cv_q_m)/(2*0.05);
    
    for l = q:num_p
        disp([q,l]);
        if l == q
            m(q,l,1) =  jac(q,1)^2; 
            m(q,l,2) =  jac(q,2)^2;
            m(q,l,3) =  jac(q,3)^2;
        else
            path_l = strcat(data_path,char(parameters(l)));
            load(strcat(path_l,'/ph/derivatives_workspace.mat'),'cell_v','cell_bud','size_at_birth','volume_at_START','cell_cycle_growth');
            if strcmp(which_dist,'full')
                size_l_p = cell_v(:,length(cell_v(1,:))) + cell_bud(:,length(cell_v(1,:)));
            elseif strcmp(which_dist,'birth')
                size_l_p = size_at_birth;
            elseif strcmp(which_dist,'start')
                size_l_p = volume_at_START;
            elseif strcmp(which_dist,'cyto')
                size_l_p = size_at_birth + cell_cycle_growth;
            end

            mean_l_p = mean(size_l_p);%(mean(size_l_p) - mean_wt)^2;
            std_l_p = std(size_l_p);%(std(size_l_p) - std_wt)^2;
            cv_l_p = std_l_p/mean_l_p;%(std_l_p/mean_l_p - cv_wt)^2;
            
            load(strcat(path_l,'/mh/derivatives_workspace.mat'),'cell_v','cell_bud','size_at_birth','volume_at_START','cell_cycle_growth');
            if strcmp(which_dist,'full')
                size_l_m = cell_v(:,length(cell_v(1,:))) + cell_bud(:,length(cell_v(1,:)));
            elseif strcmp(which_dist,'birth')
                size_l_m = size_at_birth;
            elseif strcmp(which_dist,'start')
                size_l_m = volume_at_START;
            elseif strcmp(which_dist,'cyto')
                size_l_m = size_at_birth + cell_cycle_growth;
            end

            mean_l_m = mean(size_l_m);%(mean(size_l_m) - mean_wt)^2;
            std_l_m = std(size_l_m);%(std(size_l_m) - std_wt)^2;
            cv_l_m =  std_l_m/mean_l_m;%(std_l_m/mean_l_m - cv_wt)^2;

            m(q,l,1) = (mean_q_p - mean_q_m)*(mean_l_p - mean_l_m)/((2*0.05)^2); %derivative wrt log(param)
            m(q,l,2) = (std_q_p - std_q_m)*(std_l_p - std_l_m)/((2*0.05)^2);
            m(q,l,3) = (cv_q_p - cv_q_m)*(cv_l_p - cv_l_m)/((2*0.05)^2);
            
            m(l,q,1) = m(q,l,1); %Symmetrical matrix
            m(l,q,2) = m(q,l,2);
            m(l,q,3) = m(q,l,3);
        end         
    end
end

eig_val_1 = eig(m(:,:,1));
eig_val_2 = eig(m(:,:,2));
eig_val_3 = eig(m(:,:,3));
[eig_vec_1, D] = eig(m(:,:,1));
[eig_vec_2, D] = eig(m(:,:,2));
[eig_vec_3, D] = eig(m(:,:,3));

save(strcat(save_path,'FIM_',which_dist,'_dist.mat'),'jac','m','eig_val_1','eig_val_2','eig_val_3','eig_vec_1','eig_vec_2','eig_vec_3')
%{
row_labels = {'START slope(M)','START offset','Post-START G1 (daughter) slope (M)','Post-START G1 (daughter) offset','Division (daughter) slope (M_{bud})',...
    'Division (daughter) slope (M_{START})','Division (daughter) offset','Constant growth rate','Post-START G1 (mother) slope (M)','Post-START G1 (mother) offset',...
    'Division (mother) slope (M_{bud})','Division (mother) slope (M_{START})','Division (mother) offset','Growth rate (mother) 1','Growth rate (mother) 2',...
    'Bud mass correction (daughter)','Bud mass correction (mother)','P1 synthesis rate G2','P2 synthesis G2','S/G2/M timer slope (M_{START})',...
    'S/G2/M timer offset','S/G2/M timer (mother) slope (M_{START})','S/G2/M timer (mother) offset'};     

row_labels = {'START slope(M)','START offset','Post-START G1 (daughter) slope (M)','Post-START G1 (daughter) offset','Division slope (M_{bud})',...
    'Division slope (M_{START})','Division offset','Constant growth rate','Post-START G1 (mother) slope (M)','Post-START G1 (mother) offset',...
    'Growth rate (mother) 1','Growth rate (mother) 2',...
    'Bud mass correction (daughter)','Bud mass correction (mother)','P1 synthesis rate G2','P2 synthesis G2','S/G2/M timer slope (M_{START})',...
    'S/G2/M timer offset','S/G2/M timer (mother) slope (M_{START})','S/G2/M timer (mother) offset'};

row_labels = {'START slope(M)','START offset','Post-START G1 (daughter) slope (M)','Post-START G1 (daughter) offset',...
    'Division offset','Constant growth rate','Post-START G1 (mother) slope (M)','Post-START G1 (mother) offset',...
    'Growth rate (mother) 1','Growth rate (mother) 2',...
    'Bud mass correction (daughter)','Bud mass correction (mother)','P1 synthesis rate G2','P2 synthesis G2','S/G2/M timer slope (M_{START})',...
    'S/G2/M timer offset','S/G2/M timer (mother) slope (M_{START})','S/G2/M timer (mother) offset'};
row_labels = {'START slope(M)','START offset','Post-START G1 (daughter) slope (M)','Post-START G1 (daughter) offset',...
    'Division offset','Constant growth rate','Post-START G1 (mother) slope (M)','Post-START G1 (mother) offset',...
    'Growth rate (mother) 1','Growth rate (mother) 2',...
    'Bud mass correction (daughter)','Bud mass correction (mother)','P1 synthesis rate G2','P2 synthesis G2','S/G2/M timer slope (M_{START})',...
    'S/G2/M timer offset','S/G2/M timer (mother) slope (M_{START})','S/G2/M timer (mother) offset','Growth function elbow 1','Growth function elbow 2'};

%}
row_labels = {'START slope(M)','START offset','Post-START G1 (daughter) slope (M)','Post-START G1 (daughter) offset','Division (daughter) slope (M_{bud})',...
    'Division (daughter) slope (M_{START})',['Division (daughter) offset' ...
    ''],'Constant growth rate','Post-START G1 (mother) slope (M)','Post-START G1 (mother) offset',...
    'Division (mother) slope (M_{bud})','Division (mother) slope (M_{START})','Division (mother) offset','Growth rate (mother) 1','Growth rate (mother) 2',...
    'Bud mass correction (daughter)','S/G2/M timer (daughter) slope (M_{START})',...
    'S/G2/M timer (daughter) offset','S/G2/M timer (mother) slope (M_{START})','S/G2/M timer (mother) offset','Bud mass correction (mother)'};

%% 
load(strcat(save_path,'FIM_',which_dist,'_dist.mat'),'jac','m','eig_val_1','eig_val_2','eig_val_3','eig_vec_1','eig_vec_2','eig_vec_3')

[sorted_jac_1, sorted_jac_idx_1] = sort(abs(jac(:,1)),'descend');
sorted_jac_labels_1 = row_labels(sorted_jac_idx_1);
figure('Name','Jacobian')
plot(abs(sorted_jac_1))
set(gca,'XTick',1:length(jac(:,1)))
set(gca,'XTickLabel',sorted_jac_labels_1)
set(gca,'XTickLabelRotation',45)
set(gca,'YScale','log')
xlim([1 num_p])
ylabel(texlabel('|df/d theta|'))
title(strcat('Jacobian of the MEAN of the ',which_dist,' size distributions'))
%saveas(gcf,strcat(save_path,'jacobian_mean_',which_dist,'.png'))
close

[sorted_jac_2, sorted_jac_idx_2] = sort(abs(jac(:,2)),'descend');
sorted_jac_labels_2 = row_labels(sorted_jac_idx_2);
figure('Name','Jacobian')
plot(abs(sorted_jac_2))
set(gca,'XTick',1:length(jac(:,2)))
set(gca,'XTickLabel',sorted_jac_labels_2)
set(gca,'XTickLabelRotation',45)
set(gca,'YScale','log')
xlim([1 num_p])
ylabel(texlabel('|df/d theta|'))
title(strcat('Jacobian of the STD of the ',which_dist,' size distributions'))
%saveas(gcf,strcat(save_path,'jacobian_std_',which_dist,'.png'))
close

[sorted_jac_3, sorted_jac_idx_3] = sort(abs(jac(:,3)),'descend');
sorted_jac_labels_3 = row_labels(sorted_jac_idx_3);
figure('Name','Jacobian')
plot(abs(sorted_jac_3))
set(gca,'XTick',1:length(jac(:,3)))
set(gca,'XTickLabel',sorted_jac_labels_3)
set(gca,'XTickLabelRotation',45)
set(gca,'YScale','log')
xlim([1 num_p])
ylabel(texlabel('|df/d theta|'))
title(strcat('Jacobian of the CV of the ',which_dist,' size distributions'))
%saveas(gcf,strcat(save_path,'jacobian_cv_',which_dist,'.png'))
close

eig_vec_1(abs(eig_vec_1)<0.05)=0;
[sorted_eig_val_1,sorted_eig_val_idx_1] = sort(abs(eig_val_1),'descend');
sorted_eig_vec_1 = eig_vec_1(:,sorted_eig_val_idx_1); 
cmap = flipud(brewermap(11,'RdBu'));
newCmap = imresize(cmap, [64, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);
heatmap = HeatMap(sorted_eig_vec_1,'RowLabels',row_labels,'ColumnLabels',...
    eig_val_1(sorted_eig_val_idx_1)/(max(abs(eig_val_1))),'ColorMap',newCmap,'Annotate',1);
ax = heatmap.plot;
set(gcf,'Position',[100 100 900 800]);
%set(ax,'FontSize',4)
set(ax,'XTickLabelRotation',45)
caxis(ax,[-1 1])
colorbar(ax)
addTitle(heatmap,strcat('Eigenvectors of the FIM of the MEAN of the ',which_dist,' size distributions'));
%saveas(gcf,strcat(save_path,'eigenvectors_mean_',which_dist,'.png'))
close

eig_vec_2(abs(eig_vec_2)<0.05)=0;
[sorted_eig_val_2,sorted_eig_val_idx_2] = sort(abs(eig_val_2),'descend');
sorted_eig_vec_2 = eig_vec_2(:,sorted_eig_val_idx_2); 
cmap = flipud(brewermap(11,'RdBu'));
newCmap = imresize(cmap, [64, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);
heatmap = HeatMap(sorted_eig_vec_2,'RowLabels',row_labels,'ColumnLabels',...
    eig_val_2(sorted_eig_val_idx_2)/(max(abs(eig_val_2))),'ColorMap',newCmap,'Annotate',1);
ax = heatmap.plot;
set(gcf,'Position',[100 100 900 800]);
%set(ax,'FontSize',4)
set(ax,'XTickLabelRotation',45)
caxis(ax,[-1 1])
colorbar(ax)
addTitle(heatmap,strcat('Eigenvectors of the FIM of the STD of the ',which_dist,' size distributions'));
%saveas(gcf,strcat(save_path,'eigenvectors_std_',which_dist,'.png'))
close

eig_vec_3(abs(eig_vec_3)<0.05)=0;
[sorted_eig_val_3,sorted_eig_val_idx_3] = sort(abs(eig_val_3),'descend');
sorted_eig_vec_3 = eig_vec_3(:,sorted_eig_val_idx_3); 
cmap = flipud(brewermap(11,'RdBu'));
newCmap = imresize(cmap, [64, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);
heatmap = HeatMap(sorted_eig_vec_3,'RowLabels',row_labels,'ColumnLabels',...
    eig_val_3(sorted_eig_val_idx_3)/(max(abs(eig_val_3))),'ColorMap',newCmap,'Annotate',1);
ax = heatmap.plot;
set(gcf,'Position',[100 100 900 800]);
%set(ax,'FontSize',4)
set(ax,'XTickLabelRotation',45)
caxis(ax,[-1 1])
colorbar(ax)
addTitle(heatmap,strcat('Eigenvectors of the FIM of the CV of the ',which_dist,' size distributions'));
%saveas(gcf,strcat(save_path,'eigenvectors_cv_',which_dist,'.png'))

