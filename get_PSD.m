
clear all; close all; clc

%%
eeglab
cd BCICIV_2a_gdf
% Note: To load gdf files you must download the EEGLAB plugin
[s, h] = sload('A01T.gdf', 0, 'OVERFLOWDETECTION:OFF');
[s_test, h_test] = sload('A01E.gdf', 0, 'OVERFLOWDETECTION:OFF');
%% Match indices

% Train Session
start_trial = find(h.EVENT.TYP == 768);
LH = find(h.EVENT.TYP == 769);
RH = find(h.EVENT.TYP == 770);
FT = find(h.EVENT.TYP == 771);
TG = find(h.EVENT.TYP == 772);

% Translate indices to position in dataframe
start_trial_t = h.EVENT.POS(start_trial);
LH_t = h.EVENT.POS(LH);
RH_t = h.EVENT.POS(RH);
FT_t = h.EVENT.POS(FT);
TG_t = h.EVENT.POS(TG);


% Test Session
start_trial_T = find(h_test.EVENT.TYP == 768);
LH_T = find(h_test.EVENT.TYP == 769);
RH_T = find(h_test.EVENT.TYP == 770);
FT_T = find(h_test.EVENT.TYP == 771);
TG_T = find(h_test.EVENT.TYP == 772);

% Translate indices to position in dataframe
start_trial_tT = h_test.EVENT.POS(start_trial_T);
LH_tT = h_test.EVENT.POS(LH_T);
RH_tT = h_test.EVENT.POS(RH_T);
FT_tT = h_test.EVENT.POS(FT_T);
TG_tT = h_test.EVENT.POS(TG_T);

%% Get Power Spectral Densities

LH_ch = cell(1);
RH_ch = cell(1);
FT_ch = cell(1);
TG_ch = cell(1);

% Get PSD after time delay
fs = 250; % 250 Hz sampling rate

% 1 second delay
t_start = round(fs * 1);
% 3 second delay
t_end = round(fs * 3);

for i = 1:length(LH_t)
    
    % Find the index immediately after said trial
    end_LH = min(find(start_trial_t > LH_t(i)));
    end_RH = min(find(start_trial_t > RH_t(i)));
    end_FT = min(find(start_trial_t > FT_t(i)));
    end_TG = min(find(start_trial_t > TG_t(i)));
    
    % For each trial, expand cell
    LH_ch{end+1} = s(LH_t(i)+t_start:LH_t(i)+t_end,:);
    RH_ch{end+1} = s(RH_t(i)+t_start:RH_t(i)+t_end,:);
    FT_ch{end+1} = s(FT_t(i)+t_start:FT_t(i)+t_end,:);
    TG_ch{end+1} = s(TG_t(i)+t_start:TG_t(i)+t_end,:);

end

LH_ch(1) = [];
RH_ch(1) = [];
FT_ch(1) = [];
TG_ch(1) = [];

LH_p = cell(2,1);
RH_p = cell(2,1);
FT_p = cell(2,1);
TG_p = cell(2,1);



% Compute PSD

for i = 1:length(LH_ch)
    
    x = LH_ch{i};
    [Pxx,F] = periodogram(x,[],length(x),fs);
    PSD = 10*log10(Pxx);
    LH_p{i,1} = PSD;
    LH_p{i,2} = F;
    
    x = RH_ch{i};
    [Pxx,F] = periodogram(x,[],length(x),fs);
    PSD = 10*log10(Pxx);
    RH_p{i,1} = PSD;
    RH_p{i,2} = F;
    
    x = FT_ch{i};
    [Pxx,F] = periodogram(x,[],length(x),fs);
    PSD = 10*log10(Pxx);
    FT_p{i,1} = PSD;
    FT_p{i,2} = F;
    
    x = TG_ch{i};
    [Pxx,F] = periodogram(x,[],length(x),fs);
    PSD = 10*log10(Pxx);
    TG_p{i,1} = PSD;
    TG_p{i,2} = F;
end

%% Test Trials PSD

testX_ch = cell(1);

% Get PSD after time delay
fs = 250; % 250 Hz sampling rate
% 1 second delay
t_start = round(fs * 1);
% 3 second delay
t_end = round(fs * 3);

for i = 1:length(start_trial_tT)
    
    if i == 288
        % For each trial, expand cell
        testX_ch{end+1} = s_test(start_trial_tT(i)+t_start:start_trial_tT(i)+t_end,:);
    else
        % For each trial, expand cell
        testX_ch{end+1} = s_test(start_trial_tT(i)+t_start:start_trial_tT(i)+t_end,:);
    end
    
end

testX_ch(1) = [];
testX_p = cell(2,1);

% Compute PSD

for i = 1:length(testX_ch)
    
    x = testX_ch{i};
    [Pxx,F] = periodogram(x,[],length(x),fs);
    PSD = 10*log10(Pxx);
    testX_p{i,1} = PSD;
    testX_p{i,2} = F;
    
end
%% Plot PSD of Subject 1's Different Images

figure
hold on 
grid on
shadedErrorBar(LH_p{1,2},LH_p{1,1}',{@median,@std},'lineProps','-r');
shadedErrorBar(RH_p{1,2},RH_p{1,1}',{@median,@std},'lineProps','-b');
shadedErrorBar(FT_p{1,2},FT_p{1,1}',{@median,@std},'lineProps','-c');
shadedErrorBar(TG_p{1,2},TG_p{1,1}',{@median,@std},'lineProps','-g');
title("Power Spectral Densities at the Beginning of Experiment")
xlabel("Frequency (Hz)")
ylabel("Power Spectral Density (W/Hz)")
legend('Left Hand', 'Right Hand', 'Foot', 'Tongue')
hold off

%% PSD using Somatosensory Cortex

figure
hold on 
grid on
shadedErrorBar(LH_p{1,2},LH_p{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-r');
shadedErrorBar(RH_p{1,2},RH_p{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-b');
shadedErrorBar(FT_p{1,2},FT_p{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-c');
shadedErrorBar(TG_p{1,2},TG_p{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-g');
title("Power Spectral Densities using Somatosensory Cortex")
xlabel("Frequency (Hz)")
ylabel("Power Spectral Density (W/Hz)")
legend('Left Hand', 'Right Hand', 'Foot', 'Tongue')
hold off

%% Bin Frequencies

LH_bin = cell(1);
RH_bin = cell(1);
FT_bin = cell(1);
TG_bin = cell(1);
testX_bin = cell(1);

for i = 1:length(LH_p)
    
    hold_LH = zeros(1,30-8+1);
    hold_RH = zeros(1,30-8+1);
    hold_FT = zeros(1,30-8+1);
    hold_TG = zeros(1,30-8+1);
    
    for j = 8:30
        
        % Find frequencies from 8-30
        freq_low_LH = min(find(LH_p{i,2} > j));
        freq_high_LH = min(find(LH_p{i,2} > j+1))-1;
        freq_low_RH = min(find(RH_p{i,2} > j));
        freq_high_RH = min(find(RH_p{i,2} > j+1))-1;
        freq_low_FT = min(find(FT_p{i,2} > j));
        freq_high_FT = min(find(FT_p{i,2} > j+1))-1;
        freq_low_TG = min(find(TG_p{i,2} > j));
        freq_high_TG = min(find(TG_p{i,2} > j+1))-1;
        
        % Bin frequency
        power_LH = mean(LH_p{i,1}(freq_low_LH:freq_high_LH,:));
        power_RH = mean(RH_p{i,1}(freq_low_RH:freq_high_RH,:));
        power_FT = mean(FT_p{i,1}(freq_low_FT:freq_high_FT,:));
        power_TG = mean(TG_p{i,1}(freq_low_TG:freq_high_TG,:));
        
        % Append frequency
        dim = 25;
        hold_LH(j-7,1:dim) = power_LH;
        hold_LH(j-7,dim+1) = j;
        hold_RH(j-7,1:dim) = power_RH;
        hold_RH(j-7,dim+1) = j;
        hold_FT(j-7,1:dim) = power_FT;
        hold_FT(j-7,dim+1) = j;
        hold_TG(j-7,1:dim) = power_TG;
        hold_TG(j-7,dim+1) = j;
        
    end
    
    LH_bin{i,1} = hold_LH(:,1:dim);
    LH_bin{i,2} = hold_LH(:,dim+1);
    RH_bin{i,1} = hold_RH(:,1:dim);
    RH_bin{i,2} = hold_RH(:,dim+1);
    FT_bin{i,1} = hold_FT(:,1:dim);
    FT_bin{i,2} = hold_FT(:,dim+1);
    TG_bin{i,1} = hold_TG(:,1:dim);
    TG_bin{i,2} = hold_TG(:,dim+1);
    
end

for i = 1:length(testX_p)
    
    hold_test = zeros(1,30-8+1);
    
    for j = 8:30
        
        % Find frequencies from 8-30
        freq_low_test = min(find(testX_p{i,2} > j));
        freq_high_test = min(find(testX_p{i,2} > j+1))-1;
        
        % Bin frequency
        power_test = mean(testX_p{i,1}(freq_low_test:freq_high_test,:));
        
        % Append frequency
        dim = 25;
        hold_test(j-7,1:dim) = power_test;
        hold_test(j-7,dim+1) = j;
        
    end
    
    testX_bin{i,1} = hold_test(:,1:dim);
    testX_bin{i,2} = hold_test(:,dim+1);
    
end
%% PSD using Somatosensory Cortex and 8-30 Hz With Binning

figure
hold on 
grid on
shadedErrorBar(LH_bin{1,2},LH_bin{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-r');
shadedErrorBar(RH_bin{1,2},RH_bin{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-b');
shadedErrorBar(FT_bin{1,2},FT_bin{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-c');
shadedErrorBar(TG_bin{1,2},TG_bin{1,1}(:,[8,10,12])',{@median,@std},'lineProps','-g');
title("Power Spectral Densities using Somatosensory Cortex with Binning")
xlabel("Frequency (Hz)")
ylabel("Power Spectral Density (W/Hz)")
legend('Left Hand', 'Right Hand', 'Foot', 'Tongue')
xlim([8 30])
hold off

%% Export train and test csvs

train_X = zeros(1,length([LH_bin{1,1}(:,1)', LH_bin{1,1}(:,3)', LH_bin{1,1}(:,5)']));
train_Y = zeros(1,1);

for i = 1:length(LH_p)
    
    train_X(end+1,:) = [LH_bin{i,1}(:,1)', LH_bin{i,1}(:,3)', LH_bin{i,1}(:,5)']; 
    train_Y(end+1) = 1;
    train_X(end+1,:) = [RH_bin{i,1}(:,1)', RH_bin{i,1}(:,3)', RH_bin{i,1}(:,5)']; 
    train_Y(end+1) = 2;
    train_X(end+1,:) = [FT_bin{i,1}(:,1)', FT_bin{i,1}(:,3)', FT_bin{i,1}(:,5)']; 
    train_Y(end+1) = 3;
    train_X(end+1,:) = [TG_bin{i,1}(:,1)', TG_bin{i,1}(:,3)', TG_bin{i,1}(:,5)'];
    train_Y(end+1) = 4;

end
train_X(1,:) = [];
train_Y(1) = [];

test_X = zeros(1,length([testX_bin{1,1}(:,1)', testX_bin{1,1}(:,3)', testX_bin{1,1}(:,5)']));
test_Y = zeros(1,1);

for i = 1:length(testX_p)
    
    test_X(end+1,:) = [testX_bin{i,1}(:,1)', testX_bin{i,1}(:,3)', testX_bin{i,1}(:,5)']; 

end
test_X(1,:) = [];

test_Y = load('A01E.mat');
test_Y = test_Y.classlabel';

writematrix(train_X,'../train_X.csv') 
writematrix(train_Y,'../train_Y.csv')
writematrix(test_X,'../test_X.csv') 
writematrix(test_Y,'../test_Y.csv') 
cd ..