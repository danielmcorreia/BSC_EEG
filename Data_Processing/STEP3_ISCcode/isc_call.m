% Adaptation of Dmochowski & Parra script for ISC with EEG data.
% DCorreia - Jan2022

addpath('npy-matlab')
npy_rejected_trials_drug1 = readNPY('/Volumes/LaCie Armaz/FINAL_PREPROCESSING_PAPER/PREPROCESSED_out/keep_track_rej_drug1.npy');
npy_rejected_trials_drug2 = readNPY('/Volumes/LaCie Armaz/FINAL_PREPROCESSING_PAPER/PREPROCESSED_out/keep_track_rej_drug2.npy');

for video=1:20     
   
    fprintf("\n Running ISC for Video%d", video)
    
    path_import = '/Volumes/LaCie Armaz/FINAL_PREPROCESSING_PAPER/ISC_STRUCT_out/';
    path_save = '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/';

    % some ISC processing parameters
    gamma = 0.5; % shrinkage parameter; smaller gamma for less regularization
    fs=250; % in Hz
    permutations=250; %number of permutations to run 250
    Nsec  = 1;  % time-window (in seconds) over which to compute time-resolved ISC
    overlap = 0.8*Nsec; %(1-(1/30))*Nsec; --> to be same as frame rate
    Ncomps = 3;

    stabilize_sec = 1; %number of seconds the signal takes to stabilize (exclude from ISC) %update 11_12_22
    
    % Load the data
    video_num = video;
    video = strcat("VIDEO",num2str(video));
    load(strcat(path_import, 'DRUG1/', video, '.mat'),'X');
    X1=X((stabilize_sec * fs + 1):end,:,:); %excludes 1st sec (stabilize sec) of array) %update 11_12_22
    load(strcat(path_import, 'DRUG2/', video, '.mat'),'X');
    X2=X((stabilize_sec * fs + 1):end,:,:); %excludes 1st sec (stabilize sec) of array) %update 11_12_22
    
    %update 1902 - keeping track of rejected/kept individuals

    % get arrays of equal size of individuals for group1 and group2
    min_individuals = min(size(X1,3), size(X2,3)); %gets minimum size of both groups to equalize number of individuals
    
    % gets array with original length
    X1_choose = 1:size(npy_rejected_trials_drug1,2); X2_choose = 1:size(npy_rejected_trials_drug2,2);

    % adds status of current rejections in second line
    X1_choose(2,:) = npy_rejected_trials_drug1(video_num,:); X2_choose(2,:) = npy_rejected_trials_drug2(video_num,:);
    
    % set the elements corresponding to the second row equal to 1 to their cumulative sum, the rest NAN
    X1_choose(3, X1_choose(2, :) == 1) = 1; X1_choose(3, :) = cumsum(X1_choose(3, :) == 1); X1_choose(3, X1_choose(2, :) == 0) = NaN;
    X2_choose(3, X2_choose(2, :) == 1) = 1; X2_choose(3, :) = cumsum(X2_choose(3, :) == 1); X2_choose(3, X2_choose(2, :) == 0) = NaN;
    
    % gets indexes to keep from available individuals
    idx_kept_X1 = datasample(rmmissing(X1_choose(3,:)),min_individuals,'Replace', false);
    idx_kept_X2 = datasample(rmmissing(X2_choose(3,:)),min_individuals,'Replace', false);
    
    % gets data to keep
    X1 = X1(:,:,idx_kept_X1); X2 = X2(:,:,idx_kept_X2);
    
    % Replace selected individuals by 2 to keep track
    X1_choose(2, ismember(X1_choose(3, :), idx_kept_X1)) = 2;
    X2_choose(2, ismember(X2_choose(3, :), idx_kept_X2)) = 2;

    npy_rejected_trials_drug1(video_num, :) = X1_choose(2,:); npy_rejected_trials_drug2(video_num, :) = X2_choose(2,:);
    
    %

    % ISC DRUG1
    [ISC_drug1, W_drug1, A_drug1, ISC_persubject_drug1, ISC_persubject_sum_drug1, ISC_persecond_drug1] = isc(X1, gamma, fs, Nsec, overlap);

    % ISC DRUG2
    [ISC_drug2, W_drug2, A_drug2, ISC_persubject_drug2, ISC_persubject_sum_drug2, ISC_persecond_drug2] = isc(X2, gamma, fs, Nsec, overlap);

    % calculates difference between timeseries ISC_DRUG1 - ISC_DRUG2
    difference_ISC_persecond = ISC_persecond_drug1 - ISC_persecond_drug2;
    
    save(strcat(path_save, 'DRUG1/', video, '_spatial_filter_eigval'), 'ISC_drug1');
    save(strcat(path_save, 'DRUG2/', video, '_spatial_filter_eigval'), 'ISC_drug2');

    save(strcat(path_save, 'DRUG1/', video, '_eigenvectors'), 'W_drug1');
    save(strcat(path_save, 'DRUG2/', video, '_eigenvectors'), 'W_drug2');

    save(strcat(path_save, 'DRUG1/', video, '_scalp_projection'), 'A_drug1');
    save(strcat(path_save, 'DRUG2/', video, '_scalp_projection'), 'A_drug2');

    save(strcat(path_save, 'DRUG1/', video, '_ISC_persubject'), 'ISC_persubject_drug1');
    save(strcat(path_save, 'DRUG2/', video, '_ISC_persubject'), 'ISC_persubject_drug2');

    save(strcat(path_save, 'DRUG1/', video, '_ISC_persubject_sum'), 'ISC_persubject_sum_drug1');
    save(strcat(path_save, 'DRUG2/', video, '_ISC_persubject_sum'), 'ISC_persubject_sum_drug2');

    save(strcat(path_save, 'DRUG1/', video, '_ISC_persecond'), 'ISC_persecond_drug1');
    save(strcat(path_save, 'DRUG2/', video, '_ISC_persecond'), 'ISC_persecond_drug2');
   
    % Permutation Tests
    % Check for statistical significance
    
    fprintf("\nPermuting times %d", permutations)
    
    chance_val=[];

    for iter=1:permutations
        
        % randomize data
        Xr1=phase_randomized(X1);
        Xr2=phase_randomized(X2);
        
        % run ISC code for randomized data
        [rISC_drug1, rW_drug1, rA_drug1, rISC_persubject_drug1, rISC_persubject_sum_drug1, rISC_persecond_drug1] = isc(Xr1, gamma, fs, Nsec, overlap);
        [rISC_drug2, rW_drug2, rA_drug2, rISC_persubject_drug2, rISC_persubject_sum_drug2, rISC_persecond_drug2] = isc(Xr2, gamma, fs, Nsec, overlap);
        
        chance_val_eig1(:,:,iter)=rISC_drug1';
        chance_val_eig2(:,:,iter)=rISC_drug2';

        chance_val_persubject1(:,:,iter)=rISC_persubject_drug1';
        chance_val_persubject2(:,:,iter)=rISC_persubject_drug2';

        chance_val_persubject_sum1(:,:,iter)=rISC_persubject_sum_drug1';
        chance_val_persubject_sum2(:,:,iter)=rISC_persubject_sum_drug2';

        chance_val_perwindow1(:,:,iter)=rISC_persecond_drug1';
        chance_val_perwindow2(:,:,iter)=rISC_persecond_drug2';

        % calculates difference between timeseries ISC_DRUG1 - ISC_DRUG2
        chance_vals_differences(:,:,iter) = chance_val_perwindow1(:,:,iter) - chance_val_perwindow2(:,:,iter);
        
        fprintf("\nfinished_iter%d", iter)

    end
    
    save(strcat(path_save, 'DRUG1/', video, '_chance_spatial_filter_eigval'), 'chance_val_eig1');
    save(strcat(path_save, 'DRUG2/', video, '_chance_spatial_filter_eigval'), 'chance_val_eig2');

    save(strcat(path_save, 'DRUG1/', video, '_chance_val_persubject'), 'chance_val_persubject1');
    save(strcat(path_save, 'DRUG2/', video, '_chance_val_persubject'), 'chance_val_persubject2');

    save(strcat(path_save, 'DRUG1/', video, '_chance_val_persubject_sum'), 'chance_val_persubject_sum1');
    save(strcat(path_save, 'DRUG2/', video, '_chance_val_persubject_sum'), 'chance_val_persubject_sum2');

    save(strcat(path_save, 'DRUG1/', video, '_chance_val_perwindow'), 'chance_val_perwindow1');
    save(strcat(path_save, 'DRUG2/', video, '_chance_val_perwindow'), 'chance_val_perwindow2');
    
    fprintf('\nGetting the p-values...')
  
    % for difference
    [pvals, pval_adjust]=get_pvals(chance_vals_differences, difference_ISC_persecond);
    save(strcat(path_save, 'stats/pvals_difference/', video, '_pvals_adjust_difference'), 'pval_adjust');
    save(strcat(path_save, 'stats/pvals_difference/', video, '_pvals_difference'), 'pvals');

    % for Drug1 significance persecond
    [pvals, pval_adjust]=get_pvals(chance_val_perwindow1, ISC_persecond_drug1);
    save(strcat(path_save, 'DRUG1/', video, '_pvals_adjust_signific_persec'), 'pval_adjust');
    save(strcat(path_save, 'DRUG1/', video, '_pvals_signific_persec'), 'pvals');

    % for Drug2 significance persecond
    [pvals, pval_adjust]=get_pvals(chance_val_perwindow2, ISC_persecond_drug2);
    save(strcat(path_save, 'DRUG2/', video, '_pvals_adjust_signific_persec'), 'pval_adjust');
    save(strcat(path_save, 'DRUG2/', video, '_pvals_signific_persec'), 'pvals');

    % for Drug1 significance persubject
    [pvals, pval_adjust]=get_pvals(chance_val_persubject1, ISC_persubject_drug1);
    save(strcat(path_save, 'DRUG1/', video, '_pvals_adjust_signific_persubject'), 'pval_adjust');
    save(strcat(path_save, 'DRUG1/', video, '_pvals_signific_persubject'), 'pvals');

    % for Drug2 significance persubject
    [pvals, pval_adjust]=get_pvals(chance_val_persubject2, ISC_persubject_drug2);
    save(strcat(path_save, 'DRUG2/', video, '_pvals_adjust_signific_persubject'), 'pval_adjust');
    save(strcat(path_save, 'DRUG2/', video, '_pvals_signific_persubject'), 'pvals');
    
    % for Drug1 significance persubject_sum
    [pvals, pval_adjust]=get_pvals(repmat(chance_val_persubject_sum1,1,3), repmat(ISC_persubject_sum_drug1,3,1));
    pvals = pvals(:,1); pval_adjust = pval_adjust(:,1);
    save(strcat(path_save, 'DRUG1/', video, '_pvals_adjust_signific_persubject_sum'), 'pval_adjust');
    save(strcat(path_save, 'DRUG1/', video, '_pvals_signific_persubject_sum'), 'pvals');

    % for Drug2 significance persubject_sum
    [pvals, pval_adjust]=get_pvals(repmat(chance_val_persubject_sum2,1,3), repmat(ISC_persubject_sum_drug2,3,1));
    pvals = pvals(:,1); pval_adjust = pval_adjust(:,1);
    save(strcat(path_save, 'DRUG2/', video, '_pvals_adjust_signific_persubject_sum'), 'pval_adjust');
    save(strcat(path_save, 'DRUG2/', video, '_pvals_signific_persubject_sum'), 'pvals');

    % for Drug1 significance spatial_filter_eigval
    [pvals, pval_adjust]=get_pvals(chance_val_eig1, ISC_drug1);
    save(strcat(path_save, 'DRUG1/', video, '_pvals_adjust_signific_eigen'), 'pval_adjust');
    save(strcat(path_save, 'DRUG1/', video, '_pvals_signific_eigen'), 'pvals');

    % for Drug2 significance spatial_filter_eigval
    [pvals, pval_adjust]=get_pvals(chance_val_eig2, ISC_drug2);
    save(strcat(path_save, 'DRUG2/', video, '_pvals_adjust_signific_eigen'), 'pval_adjust');
    save(strcat(path_save, 'DRUG2/', video, '_pvals_signific_eigen'), 'pvals');

    % clear
    clearvars -except npy_rejected_trials_drug1 npy_rejected_trials_drug2

end

writeNPY(npy_rejected_trials_drug1, '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/keep_track_rej_drug1_afterISC.npy');
writeNPY(npy_rejected_trials_drug2, '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/keep_track_rej_drug2_afterISC.npy');