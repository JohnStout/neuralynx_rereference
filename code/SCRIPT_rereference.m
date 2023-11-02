%% Common average rereference (CAR) CSCs, then save them back as .ncs files
%
% This code only works on PC and requires that you have included the
% Nlx2Mat .mex compiled code. Just ask neuralynx or retrieve the Nlx2Mat
% folder from https://github.com/GriffinLabCode/GriffinCode/tree/master/1.%20Matlab%20Pipeline/1.%20Formatting%20Data/Nlx2Mat
%
% Written by John Stout - edit on 11/2/23

% mode for rereferencing (average or median)
mode_reref = 'average'; % alternative is 'median'

% hardcode this
nlx_folder = getCurrentPath;
addpath(nlx_folder);

% Written by John Stout - 9/27/23
datafolder = input('Enter directory of data: ','s');
cd(datafolder);

% csc_names can be of type str or double
% csc_names are the names of your csc data. If 1-16, just use double.
% If csc_names are like 'TT1a', then use str
csc_names = [{'TT1a'} {'TT1b'} {'TT1c'} {'TT1d'} ...
    {'TT2a'} {'TT2b'} {'TT2c'} {'TT2d'} ...
    {'TT3a'} {'TT3b'} {'TT3c'} {'TT3d'} ...
    {'TT4a'} {'TT4b'} {'TT4c'} {'TT4d'}];

% this is when the LFPs are all named 'CSC'
%csc_names = [1:16];

%% Load in CSC data
disp('Loading in CSC channels. This could take a while...')
Samples3D = [];
for i = 1:length(csc_names)

    % define csc name in raw format
    if contains(class(csc_names),'double')
        varName = strcat('\CSC',num2str(csc_names(i)),'.ncs');
    else
        varName = strcat('\',csc_names{i},'.ncs');    
    end
    
    % define directory
    dir = strcat(datafolder,varName);

    % convert CSC
    if contains(class(csc_names),'double')
        disp(['loading CSC',num2str(csc_names(i))]);
    else
        disp(['loading ', csc_names{i}])
    end
    
    % Neuralynx function
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
        Samples, Header] = Nlx2MatCSC(dir, [1 1 1 1 1], 1, 1, []);

    % concatenate
    Samples3D(:,:,i) = Samples;
end

%% visualize briefly

% Get the center of the dataset and plot 10s
idx_mid = round(length(Timestamps)/2);
srate_crunch = round(32000/512);
ts_start = find(Timestamps==Timestamps(idx_mid-(srate_crunch*10)));
ts_end = find(Timestamps==Timestamps(idx_mid+(srate_crunch*10)));

% colors
tt_counts = size(Samples3D,3)/4;
channel_counts = 4;
colors = [];
for i = 1:length(tt_counts)
    colors{i} = [round(rand) round(rand) round(rand)];
end

disp('This figure assumes that every 4 wires is a new tetrode')
tt_counter = 1; rem_ch = []; channel_counter = [];
for i = 1:4:size(Samples3D,3)
    figure('color','w');
    channel_counter{tt_counter} = i:i+3;
    data_tt = []; data_channel = [];
    data_tt = Samples3D(:,ts_start:ts_end,i:i+3); % samples by tetrode from center point
    for ploti = 1:size(data_tt,3)
        % Now lets extrac channel data
        temp_data = [];
        temp_data = data_tt(:,:,ploti); % temp variable
        data_channel(:,ploti) = temp_data(:); % vectorize and save data
    end
    % plot results is separate to troubleshoot data_channel if needed
    for ploti = 1:size(data_tt,3)
        % plot data
        subplot(tt_counts,1,ploti)
        plot(data_channel(:,ploti),'k')
        title(['Tetrode',num2str(tt_counter),', Channel',num2str(ploti)])
    end
    rem_ch{tt_counter} = str2num(input('Enter channels to exclude from CAR (e.g. enter "1,3" to remove channels 1 and 3): ','s'));
    tt_counter = tt_counter+1;
    close;
end
close all;

% convert channel counter to actual channel number (1-16, not 0-15) to
% index out the Samples3D
ch_ex = [];
for i = 1:length(rem_ch) % loop over tetrode
    ch_ex{i} = channel_counter{i}(rem_ch{i});
end
ch_ex = horzcat(ch_ex{:});

% rather than setting to NaN and eating up memory, lets just ignore
channel_idx = 1:tt_counts*4;
channel_log = logical(1:tt_counts*4);
channel_log(ch_ex)=0;
disp('Getting avg across channels. This might take a minute...')

%% mode for rereferencing
if contains(mode,'median')
    disp('Collecting the common median')
    cr = median(Samples3D(:,:,channel_log),3);
else % default is average
    disp('Collecting the common average')
    cr = mean(Samples3D(:,:,channel_log),3);
end

%% Common average rereference
clearvars -except datafolder csc_names Samples3D ca channel_idx channel_log ch_idx
disp('Performing the common avg rereferencing... may take a few moments...')
SamplesRef = zeros(size(Samples3D));
SamplesRef = cell([1 size(Samples3D,3)]);
for i = 1:size(Samples3D,3)
    SamplesRef{i} = Samples3D(:,:,i)-cr;
end
% remove tagged channels
SamplesRef{~channel_log}=[];

%% cross correlation
data2cor_og = []; data2cor_ref = [];
for i = 1:size(Samples3D,3)
    % first og data
    temp_data = []; temp_data = Samples3D(:,1:500,i);
    temp_data = temp_data(:);
    data2cor_og(:,i) = temp_data;
    
    % then rereferenced data
    temp_data = []; 
    if isempty(SamplesRef{i})
        temp_data = NaN(size(Samples3D(:,1:500,1)));
    else
        temp_data = SamplesRef{i}(:,1:500);
    end
    temp_data = temp_data(:);    
    data2cor_ref(:,i) = temp_data;    
end
R_og = corrcoef(data2cor_og);
R_ref = corrcoef(data2cor_ref);
figure('color','w')
subplot(2,1,1);
    pcolor(R_og);
    colorbar;
    %shading interp
    caxis([-1 1])
    %clabel('R')
    title('OG signal')
    ylabel('Channel #')
subplot(2,1,2);
    pcolor(R_ref);
    colorbar;
    %shading interp
    caxis([-1 1])
    title('Reref signal')
    xlabel('Channel #')
    %clabel('R')
savefig('fig_crosscorr_rereference');

%% Figure showing difference between OG signal and rereferenced signal
disp('INCOMPLETE - MUST MARK DEAD/EXCLUDED CHANNELS')
figure('color','w'); counter = 1; fs = 32000;
for i = 1:size(data2cor_og,2)
    subplot(size(data2cor_og,2),1,counter); box off;
    plot(data2cor_og(1:fs*3,i),'k'); hold on;
    plot(data2cor_ref(1:fs*3,i),'r');
    counter = counter+1;
end
legend('OG signal','Reref signal')
savefig('fig_channelComparison_rereference');

%% Save out CSC data
disp('Saving output. This could take a minute...'); counter = 1;
for i = csc_names

    % define csc name in raw format
    if contains(class(csc_names),'double')
        varName = strcat('\CSC',num2str(csc_names(i)),'.ncs');
    else
        varName = strcat('\',csc_names{i},'.ncs');    
    end
    
    % checking for data and skipping if not present
    if isempty(SamplesRef{i})
        disp(['Skipping ',varName])
        continue
    end

    % load csc
    disp(['loading csc',num2str(i)])
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
        ~, Header] = Nlx2MatCSC(strcat(datafolder,varName), [1 1 1 1 1], 1, 1, []);

    % use the loaded data to save out csc
    varSave = ['\csc_car',num2str(i),'.ncs'];
    AppendToFileFlag = 0; 
    ExportMode = 1;
    ExportModeVector = 1;
    FieldSelectionFlags = [1 1 1 1 1 1];

    % write csc file
    Samples = [];
    Samples = SamplesRef{i};
    disp(['Writing csc_ca',num2str(i),'.ncs'])
    Mat2NlxCSC(strcat(datafolder,varSave), AppendToFileFlag, ExportMode, ExportModeVector,...
         FieldSelectionFlags, Timestamps, ChannelNumbers,...
         SampleFrequencies, NumberOfValidSamples, Samples, Header);
     
    % add to counter
    counter = counter+1;

end

%{
%% troubleshooting purposes (ignore)
[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
    Samples, Header] = Nlx2MatCSC(strcat(datafolder,'/CSC1.ncs'), [1 1 1 1 1], 1, 1, []);
csc1 = Samples(:)';

[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
    Samples, Header] = Nlx2MatCSC(strcat(datafolder,'/csc_car1.ncs'), [1 1 1 1 1], 1, 1, []);
csc1_car = Samples(:)';

figure('color','w'); plot(csc1(1:10000),'b'); hold on; plot(csc1_car(1:10000),'r');
legend('OG signal','CAR signal'); box off;

fs = 32000;
csc1_filt = skaggs_filter_var(csc1,500,9000,32000);
csc1_car_filt = skaggs_filter_var(csc1_car,500,9000,32000);

figure; plot(csc1_filt(1:5000),'b');
hold on; plot(csc1_car_filt(1:5000),'r');
legend('OG signal','CAR signal'); box off;
%}
