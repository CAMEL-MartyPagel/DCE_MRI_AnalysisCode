% This code is written by Chetan B Dhakan on 08/25/2025
% This code works for DCE-MRI scanned from Bruker preclinical MRI scanners
% It populates Concentration, R1 and delta R1 over time
% It also computes Ktrans, Kep from Tofts model and RKtrans, Kep,
% Kepmuscle from LRRM model and also populates their parametric maps.
%% Clear Workspace and define pathnames

clearvars; clc; close all;
Starting_Directory=pwd;
%% Identify directory for saving results. Load anatomical and DCE pre and post images

cd(Starting_Directory);
disp('Select the directory for saving the final results ');
DIR=uigetdir; 
parent_DIR=strcat(DIR,'/');

cd(parent_DIR);
disp('Select anatomical RARE scan ');
DIR=uigetdir; 
DIR_anatomical=strcat(DIR,'/');

[anatomical_image,anatomical_offsets,TR ] = load_anatomical_image(DIR_anatomical);
dim = size(anatomical_image);

cd(parent_DIR);
disp('Select DCE pre scan ');
DIR=uigetdir; 
DIR_DCE=strcat(DIR,'/');

cd(parent_DIR);
disp('Select DCE post scan ');
DIR=uigetdir; 
DIR_DCE_post=strcat(DIR,'/');

[DCE_image, DCE_TR, flip_angle] = load_DCE_image(DIR_DCE);

[DCE_image_post, DCE_TR_post, flip_angle_post] = load_DCE_image(DIR_DCE_post);

%% Number of ROIs to analyze
number_of_ROIs = 3;  % hard coding; Change this to any number if you want to analyze more than two ROIs

% Initialize a cell array to store ROI results, here 35 and 135 are hard
% coded but change it to your DCE MRI no. of scans pre and post.

R1_ROI=zeros(number_of_ROIs,35);
R1_ROI_post = zeros(number_of_ROIs, 135);
conc = zeros(number_of_ROIs,35);
conc_post = zeros(number_of_ROIs, 135);

fh = figure;
imagesc(anatomical_image(:,:,1,1));colormap('gray');axis tight;
title (['anatomical image']);
exportgraphics(fh, "anatomicalimage.jpg", "Resolution", 300);

for n=1:3 %number_of_ROIs
    if n == 1
        %prompt = 'What is the name of the first ROI?';
        name_of_region = 'tumor';
    elseif n == 2
        %prompt = 'What is the name of the second ROI?';
        name_of_region = 'muscle';
    elseif n == 3
        %prompt = 'What is the name of the third ROI?';
        name_of_region = 'aif';
%     else
%         prompt = sprintf('What is the name of ROI #%d?', n);
    end

    
    %name_of_region = inputdlg(prompt, 'Input ROI Name', [1 50]);
        
    f1 = figure;
    imagesc(anatomical_image(:,:,1,1));colormap('gray');axis tight;
    title(['Choose ',name_of_region, ' ROI ']);
    h = drawfreehand;
    mask = createMask(h);
    close(f1);

    if n == 1
        mask_tumor = mask;
    elseif n == 2
        mask_muscle = mask;
    elseif n == 3
        mask_aif = mask;
    end

    masked_anatomical = squeeze(anatomical_image(:,:,1,:).*mask);
    DCE_slice_selected = squeeze(DCE_image(:,:,1,:).*mask);
    DCE_slice_selected_post = squeeze(DCE_image_post(:,:,1,:).*mask);
        
    T1map = zeros(dim(1),dim(2));
    M0 = zeros(dim(1),dim(2));
    c=zeros(dim(1),dim(2));
    
    for i = 1:dim(1)
        for j = 1:dim(2)
            if masked_anatomical(i,j,1)~=0
                [M0(i,j), T1map(i,j),c(i,j)] = t1fitting_VTR(transpose(squeeze(masked_anatomical(i,j,:))), TR);
                %[M0(i,j), T1map(i,j),c(i,j)] = t1fitting_VFA(transpose(squeeze(masked_anatomical(i,j,:))), TR,FAs);
            end
        end
    end
    
    % 35 timepoints were acquired prior to injection of contrast; 
    DCE_before = DCE_slice_selected(:,:,1,:); 
    DCE_after = DCE_slice_selected_post(:,:,1,:);
    % S0: pre-contrast signal intensity. 
    S0 = mean(DCE_before,3);
    S = mean(DCE_after,3);
    
    dim_DCE=size(DCE_slice_selected);
    dim_DCE_post = size(DCE_slice_selected_post);
    
    % R1: Relaxivity and R1_ROI: avergage relaxivity within ROI.
    R1=zeros(dim_DCE(1),dim_DCE(2),dim_DCE(3));
    R1_post = zeros(dim_DCE_post(1),dim_DCE_post(2),dim_DCE_post(3));
    T1_0 = 0;
    T1_0_post = 0;
    
    for k = 1:dim_DCE(3)
        for i=1:dim_DCE(1)
            for j=1:dim_DCE(2)
                if DCE_slice_selected(i,j,k)~=0 && T1map(i,j)>0 
                    E = exp(-DCE_TR/T1map(i,j));
                    B = (1-E)/(1-(cosd(flip_angle)*E));
                    A = B*(DCE_slice_selected(i,j,k)/S0(i,j));
                    C = (1-A)/(1-A*cosd(flip_angle));
                    if C>0
                        R1(i,j,k) = log(C)*(-1/DCE_TR);
                        if R1(i,j,k) >= 0  % Only consider valid R1 values
                           T1_0 = T1_0 + T1map(i,j);
                        else
                            R1(i,j,k) = NaN;  % Set invalid R1 to N
                        end
                    end
                end
            end
        end
    end

    for k = 1:dim_DCE_post(3)
        for i=1:dim_DCE_post(1)
            for j=1:dim_DCE_post(2)
                if DCE_slice_selected_post(i,j,k)~=0 && T1map(i,j)>0 
                    E = exp(-DCE_TR_post/T1map(i,j));
                    B = (1-E)/(1-(cosd(flip_angle_post)*E));
                    A = B*(DCE_slice_selected_post(i,j,k)/S(i,j));
                    C = (1-A)/(1-A*cosd(flip_angle_post));
                    if C>0
                        R1_post(i,j,k) = log(C)*(-1/DCE_TR_post);
                        if R1_post(i,j,k) >= 0  % Only consider valid R1 values
                            T1_0_post = T1_0_post + T1map(i,j);
                        else
                            R1_post(i,j,k) = NaN;  % Set invalid R1 to NaN
                        end
                    end
                end
            end
        end
    end


    for k = 1:dim_DCE(3)
        R1_ROI(n,k)=mean(nonzeros(R1(:,:,k)));
        R1_ROI_4D(:,:,k,n) = R1(:,:,k);
    end
    for k = 1:dim_DCE_post(3)
        R1_ROI_post(n,k)=mean(nonzeros(R1_post(:,:,k)));
        R1_ROI_post_4D(:,:,k,n) = R1_post(:,:,k);
    end

    %R1_ROI{n,:}=1/T1_0;
    R1_new=1000*R1_ROI;
    R1_new_post=1000*R1_ROI_post;

    R1_new_4D=1000*R1_ROI_4D;
    R1_new_post_4D=1000*R1_ROI_post_4D;

    concentration = zeros(dim_DCE(1), dim_DCE(2), dim_DCE(3));
    concentration_post = zeros(dim_DCE_post(1), dim_DCE_post(2), dim_DCE_post(3));
    r1 = 3.8/1000;

    % Calculate concentration for each time point
for k = 1:dim_DCE(3)
    for i = 1:dim_DCE(1)
        for j = 1:dim_DCE(2)
            if DCE_slice_selected(i,j,k) ~= 0 && T1map(i,j) > 0
                % Calculate concentration using the formula
                concentration(i,j,k) = (R1(i,j,k) - (1/T1_0)) / r1; % r1 = 5
            end
        end
    end
end

% Calculate post-contrast concentration
for k = 1:dim_DCE_post(3)
    for i = 1:dim_DCE_post(1)
        for j = 1:dim_DCE_post(2)
            if DCE_slice_selected_post(i,j,k) ~= 0 && T1map(i,j) > 0
                % Calculate concentration using the formula
                concentration_post(i,j,k) = (R1_post(i,j,k) - (1/T1_0_post)) / r1; % r1 = 5
            end
        end
    end
end

    for k = 1:dim_DCE(3)
        conc(n,k)=mean(nonzeros(concentration(:,:,k)));
        conc_4D(:,:,k,n) = concentration(:,:,k);
    end
    for k = 1:dim_DCE_post(3)
        conc_post(n,k)=mean(nonzeros(concentration_post(:,:,k)));
        conc_post_4D(:,:,k,n) = concentration_post(:,:,k);
    end

end

%% T1 map

% T1map = zeros(dim_DCE_post(1),dim_DCE_post(2));
% M0 = zeros(dim_DCE_post(1),dim_DCE_post(2));
% c=zeros(dim_DCE_post(1),dim_DCE_post(2));
% 
% anatomical_image_T1_masked = squeeze(anatomical_image(:,:,1,:).*mask_tumor);
% 
% disp("Starting T1mapping ...");
% cd(Starting_Directory);
% for i = 1:dim_DCE_post(1)
%     for j = 1:dim_DCE_post(2)
%         A = squeeze(anatomical_image(i,j,1,:));
%         [M0(i,j), T1map(i,j)] = t1fitting_VTR(A, TR);
%     end
% end
% 
% T1map_masked = T1map.*mask_tumor;
% T1map_masked_mean = mean(T1map_masked(mask_tumor));
% disp(T1map_masked_mean);
% 
% % Display Results
% figure;
% ax1 = axes;
% imagesc(anatomical_image(:,:,1,1))
% axis image off;
% colormap(ax1, 'gray');
% ax2 = axes;
% imagesc(ax2, T1map, 'alphadata', mask_tumor);
% axis image off;
% colormap(ax2,'jet');
% caxis([0 3000]);
% ax2.Visible = 'off';
% linkprop([ax1 ax2],'Position');
% title('T1 map');
% colorbar;
% 
% % Save Data
% savefig('T1_map.fig');
% print(gcf, '-dtiff', 'T1_map.tiff');

%% delta R1 and concentration calculation

% Calculate the baseline R1 (mean of pre-contrast frames)
R1_pre_baseline_tumor = mean(R1_new(1,:), 'omitnan');  % Average of the 35 pre-contrast frames
R1_pre_baseline_muscle = mean(R1_new(2,:), 'omitnan');
R1_pre_baseline_aif = mean(R1_new(3,:), 'omitnan');

R1_pre_baseline_tumor_4D = mean(R1_new_4D(:,:,:,1),3, 'omitnan');  % Average of the 35 pre-contrast frames
R1_pre_baseline_muscle_4D = mean(R1_new_4D(:,:,:,2),3, 'omitnan');
R1_pre_baseline_aif_4D = mean(R1_new_4D(:,:,:,3),3, 'omitnan');

% Calculate change in R1 (Delta R1)
Delta_R1_tumor = R1_new_post(1,:) - R1_pre_baseline_tumor;
Delta_R1_muscle = R1_new_post(2,:) - R1_pre_baseline_muscle;
Delta_R1_aif = R1_new_post(3,:) - R1_pre_baseline_aif;

Delta_R1_tumor_4D = R1_new_post_4D(:,:,:,1) - R1_pre_baseline_tumor_4D(:,:);
Delta_R1_muscle_4D = R1_new_post_4D(:,:,:,2) - R1_pre_baseline_muscle_4D(:,:);
Delta_R1_aif_4D = R1_new_post_4D(:,:,:,3) - R1_pre_baseline_aif_4D(:,:);

% Extract the first row from both matrices
data_new = R1_new(1, :); % 1x35
data_post = R1_new_post(1, :); % 1x135
data_new_muscle = R1_new(2, :); % 1x35
data_post_muscle = R1_new_post(2, :); % 1x135
data_new_aif = R1_new(3, :); % 1x35
data_post_aif = R1_new_post(3, :); % 1x135

data_new_4D = R1_new_4D(:,:,:,1); % 128x128x35
data_post_4D = R1_new_post_4D(:,:,:,1); % 128x128x135
data_new_muscle_4D = R1_new_4D(:,:,:,2); % 128x128x35
data_post_muscle_4D = R1_new_post_4D(:,:,:,2); % 128x128x135
data_new_aif_4D = R1_new_4D(:,:,:,3); % 128x128x35
data_post_aif_4D = R1_new_post_4D(:,:,:,3); % 128x128x135

% Combine the data
combined_data_tumor = [data_new, data_post]; % 1x170
combined_data_muscle = [data_new_muscle, data_post_muscle];
combined_data_aif = [data_new_aif, data_post_aif];

combined_data_tumor_4D = cat(3, data_new_4D, data_post_4D); % 128x128x170
combined_data_muscle_4D = cat(3, data_new_muscle_4D, data_post_muscle_4D);
combined_data_aif_4D = cat(3, data_new_aif_4D, data_post_aif_4D);

conc_tumor = conc(1,:);
conc_tumor_post = conc_post(1,:);
conc_muscle = conc(2,:);
conc_muscle_post = conc_post(2,:);
conc_aif = conc(3,:);
conc_aif_post = conc_post(3,:);

conc_tumor_4D = conc_4D(:,:,:,1);
conc_tumor_post_4D = conc_post_4D(:,:,:,1);
conc_muscle_4D = conc_4D(:,:,:,2);
conc_muscle_post_4D = conc_post_4D(:,:,:,2);
conc_aif_4D = conc_4D(:,:,:,3);
conc_aif_post_4D = conc_post_4D(:,:,:,3);

conc_total_tumor = [conc_tumor, conc_tumor_post];
conc_total_muscle = [conc_muscle, conc_muscle_post];
conc_total_aif = [conc_aif, conc_aif_post];

conc_total_tumor_4D = cat(3, conc_tumor_4D, conc_tumor_post_4D);
conc_total_muscle_4D = cat(3, conc_muscle_4D, conc_muscle_post_4D);
conc_total_aif_4D = cat(3, conc_aif_4D, conc_aif_post_4D);

% Number of frames
numFrames = 135;

% Duration of each frame in seconds
frameDuration = 10;

% Create time vector in seconds
timeVectorSeconds = (0:numFrames-1) * frameDuration;

% Convert time vector to minutes
time = timeVectorSeconds / 60;

% Number of frames
numFrames_post = 170;

% Duration of each frame in seconds
frameDuration_post = 10;

% Create time vector in seconds
timeVectorSeconds_post = (0:numFrames_post-1) * frameDuration_post;

% Convert time vector to minutes
time_post = timeVectorSeconds_post / 60;

cd(parent_DIR);

% Plot Delta R1 vs. Time
f2 = figure;
plot(time, Delta_R1_tumor, 'ro', 'LineWidth', 1.5);
hold on
plot(time, Delta_R1_muscle, 'ko', 'LineWidth', 1.5);
plot(time, Delta_R1_aif, 'mo', 'LineWidth', 1.5);
hold off
xlabel('Time (minutes)');
ylabel('\Delta R1 (s^{-1})');
title(['\Delta R1 over Time']);
xlim([0 25]);
xticks(0:5:25);
legend('tumor', 'muscle', 'aif');
exportgraphics(f2, "deltar1tumormuscle.jpg", "Resolution",300);

% Plot Concentration vs. Time
f3 = figure;
plot(time_post, conc_total_tumor, 'ro', 'LineWidth', 1.5);
hold on
plot(time_post, conc_total_muscle, 'ko', 'LineWidth', 1.5);
plot(time_post, conc_total_aif, 'mo', 'LineWidth', 1.5);
hold off
xlabel('Time (minutes)');
ylabel('Concentration (mM)');
title(['Concentration over Time']);
xlim([0 25]);
xticks(0:5:25);
legend('tumor', 'muscle', 'aif');
exportgraphics(f3, "conctumormuscle.jpg", "Resolution",300);

f4 = figure;
plot(time_post, combined_data_tumor, 'ro', 'LineWidth', 1.5);
hold on
plot(time_post, combined_data_muscle, 'ko', 'LineWidth', 1.5);
plot(time_post, combined_data_aif, 'mo', 'LineWidth', 1.5);
hold off
xlabel('Time (minutes)');
xlim([0 25]);
xticks(0:5:25);
ylabel('R1 (1/s)');
title(['R1 over Time']);
legend('tumor', 'muscle', 'aif');
exportgraphics(f4, "R1.jpg", "Resolution",300);

%% Compute DCE-MRI kinetic parameters from Tofts and LRRM model
dcemrinonneg = fitdcemri(Delta_R1_tumor,Delta_R1_muscle,time,'nonneg');
dcemrinoneg_c = fitdcemri(conc_total_tumor,conc_total_muscle,time_post,'nonneg');

dcemrilsq = fitdcemri(Delta_R1_tumor,Delta_R1_muscle,time,'lsq');
dcemrilsq_c = fitdcemri(conc_total_tumor,conc_total_muscle,time_post,'lsq');

dcemriaif = fitdcemri(conc_total_tumor, conc_total_aif, time_post);


datatable = table(dcemrinonneg, dcemrinoneg_c, dcemrilsq, dcemrilsq_c, ...
                  'VariableNames', {'DCE_MRI_NonNeg', 'DCE_MRI_NonNeg_C', ...
                                    'DCE_MRI_LSQ', 'DCE_MRI_LSQ_C'});
% Combine all variables into a table
combined_table = table(conc_total_tumor', conc_total_muscle', conc_total_aif', ...
                       combined_data_tumor', combined_data_muscle', combined_data_aif', ...
                       'VariableNames', {'Concentration_Tumor', 'Concentration_Muscle', ...
                                         'Concentration_AIF', 'R1Tumor', ...
                                         'R1Muscle', 'R1AIF'});

save('DCE_MRI_combined_table.mat', 'combined_table');
save('DCE_MRI_Results.mat', 'datatable', 'dcemriaif');

%% Display Results
figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, mask_tumor, 'alphadata', dcemrinoneg_c(1).*mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('RKtrans');
colorbar;

%%

% Initialize output arrays with zeros
dcemrinonneg_4D = zeros(dim_DCE_post(1), dim_DCE_post(2), 4);
dcemrinoneg_c_4D = zeros(dim_DCE_post(1), dim_DCE_post(2), 4);
dcemrilsq_4D = zeros(dim_DCE_post(1), dim_DCE_post(2), 4);
dcemrilsq_c_4D = zeros(dim_DCE_post(1), dim_DCE_post(2), 4);
dcemriaif_4D = zeros(dim_DCE_post(1), dim_DCE_post(2), 3);

for i = 1:dim_DCE_post(1)
    for j = 1:dim_DCE_post(2)
        if mask_tumor(i,j) > 0
            dcemrinonneg_4D(i,j,:) = fitdcemri(squeeze(Delta_R1_tumor_4D(i,j,:)),Delta_R1_muscle,time,'nonneg');
            dcemrinoneg_c_4D(i,j,:) = fitdcemri(squeeze(conc_total_tumor_4D(i,j,:)),conc_total_muscle,time_post,'nonneg');
    
            dcemrilsq_4D(i,j,:) = fitdcemri(squeeze(Delta_R1_tumor_4D(i,j,:)),Delta_R1_muscle,time,'lsq');
            dcemrilsq_c_4D(i,j,:) = fitdcemri(squeeze(conc_total_tumor_4D(i,j,:)),conc_total_muscle,time_post,'lsq');
    
            dcemriaif_4D(i,j,:) = fitdcemri(squeeze(conc_total_tumor_4D(i,j,:)), conc_total_aif, time_post);
        end
    end
end

%% Display Results
figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, dcemriaif_4D(:,:,1), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('Ktrans - AIF');
colorbar;
saveas(gcf, 'Ktrans_LTM_AIF.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, dcemriaif_4D(:,:,2), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('kep - AIF');
colorbar;
saveas(gcf, 'kep_LTM_AIF.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, dcemrinoneg_c_4D(:,:,1), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('RKtrans - LRRM - nonneg fit - Reference: Muscle');
colorbar;
saveas(gcf, 'RKtrans_LRRM_nonneg_refMuscle.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, dcemrinoneg_c_4D(:,:,2), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('kep - LRRM - nonneg fit - Reference: Muscle');
colorbar;
saveas(gcf, 'kep_LRRM_nonneg_refMuscle.tiff.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, dcemrilsq_c_4D(:,:,1), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('RKtrans - LRRM - LLSQ fit - Reference: Muscle');
colorbar;
saveas(gcf, 'RKtrans_LRRM_LLSQ_refMuscle.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, dcemrilsq_c_4D(:,:,2), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('kep - LRRM - LLSQ fit- Reference: Muscle');
colorbar;
saveas(gcf, 'kep_LRRM_LLSQ_refMuscle.tiff');

%% Timeline images
figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, conc_total_tumor_4D(:,:,10), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('Baseline');
colorbar;
caxis([0 0.15]);
saveas(gcf, 'conc_Baseline.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, conc_total_tumor_4D(:,:,36), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('6 Minutes');
colorbar;
caxis([0 0.15]);
saveas(gcf, 'conc_6minutes.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, conc_total_tumor_4D(:,:,48), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('8 Minutes');
colorbar;
caxis([0 0.15]);
saveas(gcf, 'conc_8minutes.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, conc_total_tumor_4D(:,:,60), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('10 Minutes');
colorbar;
caxis([0 0.15]);
saveas(gcf, 'conc_10minutes.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, conc_total_tumor_4D(:,:,120), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('20 Minutes');
colorbar;
caxis([0 0.15]);
saveas(gcf, 'conc_20minutes.tiff');

figure;
ax1 = axes;
imagesc(anatomical_image(:,:,1,1))
axis image off;
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, conc_total_tumor_4D(:,:,150), 'alphadata', mask_tumor);
axis image off;
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('25 Minutes');
colorbar;
caxis([0 0.15]);
saveas(gcf, 'conc_25minutes.tiff');



