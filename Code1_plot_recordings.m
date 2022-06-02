t% ############################################################
% Code 1: plot recordings from Percept device json sessions
%
% Author: Matteo Vissani @BrainModulation Lab, Boston (MA), USA
% ############################################################

%% load all info

% set paths
config_paths;

% set figure properties
prettify_plots;

% identify patients and sessions
config_data;

% get clinical scales and med changes from csv file (date field in seconds)
config_clinical_info;

% save fig?
FLAG_FIGSAVE = true;


%% extract useful variables
ybocs = info_data.clinical_info.YBOCS;
madrs = info_data.clinical_info.MADRS;
sID = info_data.clinical_info.ID;
sessions_clin = info_data.clinical_info.date;
medchanges = info_data.med_changes.change;
medchanges_date = info_data.med_changes.date;

%% check clinical scales ---- %%

figure("renderer","painters","position",[300 300 500 500],"DefaultAxesFontSize",15)
scatter(ybocs, madrs,105,"k","filled")
text(ybocs, madrs + .5 , sID)
xlabel(" YBOCS ")
ylabel("MADRS ")
box off

if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_YBOCS_vs_MADRS"));
    saveas(gcf,figname,"png")
end


%% load json files from Percept device

PATH_JSON = fullfile(PATH_DATA,info_data.PATIENT,info_data.FILES);

data_all = cell(1,info_data.n_FILES);
n_events = nan(1,info_data.n_FILES);

for sess_i = 1: info_data.n_FILES
    data_all{1,sess_i} = jsondecode(fileread(PATH_JSON{sess_i}));
    try
        n_events(sess_i) = numel([data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents.EventID]);
    catch
    end
end
[~,idx_forevent] = max(n_events);



events = datetime({data_all{1, idx_forevent}.DiagnosticData.LfpFrequencySnapshotEvents.DateTime},'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
events_ID = [data_all{1, idx_forevent}.DiagnosticData.LfpFrequencySnapshotEvents.EventID];
events_names = {data_all{1, idx_forevent}.DiagnosticData.LfpFrequencySnapshotEvents.EventName};
events_unique = unique(events_ID);

switch numel(events_unique)
    case 1
        colors_event = {[0.3010 0.7450 0.9330]};
    case 2
        colors_event = {[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880]};
    case 3
        colors_event = {[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560]};
end


lfp_left = [];
lfp_right = [];
time_left = [];
time_right = [];
stim_left = [];
stim_right = [];
session = [];


for sess_i = 1: info_data.n_FILES

    session = [session datetime(data_all{1,sess_i}.SessionEndDate,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''')];
    fields_left = fieldnames(data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left);
    for field_i = 1 : numel(fields_left)
        cdt = {data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left.(fields_left{field_i}).DateTime};
        time_left = [time_left datetime(cdt,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');];
        topick = [data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left.(fields_left{field_i}).LFP];
        lfp_left = [lfp_left topick];
        stim_left = [stim_left [data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left.(fields_left{field_i}).AmplitudeInMilliAmps]];
    end

    fields_right = fieldnames(data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Right);

    for field_i = 1 : numel(fields_right)
        lfp_right = [lfp_right [data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Right.(fields_right{field_i}).LFP]];
        stim_right = [stim_right [data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Right.(fields_right{field_i}).AmplitudeInMilliAmps]];
        cdt = {data_all{1, sess_i}.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Right.(fields_right{field_i}).DateTime};
        time_right = [time_right datetime(cdt,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''')];
    end
end

% adjust left samples order
[~,idx] = sort(datenum(time_left));
time_left_s = time_left(idx);
lfp_left_s = lfp_left(idx);
[~,idx_unique] = unique(time_left_s);
lfp_left_s = lfp_left_s(idx_unique);

stim_left_s = stim_left(idx);
stim_left_s = stim_left_s(idx_unique);
time_left_s = time_left_s(idx_unique);

% adjust right samples order
[~,idx] = sort(datenum(time_right));
time_right_s = time_right(idx);
lfp_right_s = lfp_right(idx);
[~,idx_unique] = unique(time_right_s);
lfp_right_s = lfp_right_s(idx_unique);

stim_right_s = stim_right(idx);
stim_right_s = stim_right_s(idx_unique);
time_right_s = time_right_s(idx_unique);


% correct time vector for stimulation variables because of a gap in the
% data
time_left_s_stim = [time_left_s(1:1142) datetime(2020,12,11):minutes(10):datetime(2021,1,10) time_left_s(1143:end)];
time_right_s_stim = [time_right_s(1:1142) datetime(2020,12,11):minutes(10):datetime(2021,1,10) time_right_s(1143:end)];

stim_left_s = [stim_left_s(1:1142) 4*ones(1,numel(datetime(2020,12,11):minutes(10):datetime(2021,1,10))) stim_left_s(1143:end)];
stim_right_s  = [stim_right_s(1:1142) 4*ones(1,numel(datetime(2020,12,11):minutes(10):datetime(2021,1,10)))  stim_right_s(1143:end)];



%% Plot raw data after a simple movmedian operation to get rid of brief transient artefacts

figure("renderer","painters","Position",[400   400   850 600])
t = tiledlayout(7,1);
% plot stimulation
nexttile(1,[1,1])
times =  time_left_s_stim - time_left_s(1);
times.Format = "d";
plot(times,stim_left_s,"r")
hold on
ylabel(" Stim Amp. [mA] ")
ylim([0 7])
xlim([times(1) - days(5) times(end)])
box off

% plot excessive thoughts events
nexttile(2,[1,1])
for event_i = 1: numel(events_unique)
    plot([events(events_ID == events_unique(event_i)); events(events_ID == events_unique(event_i)) ]- time_right_s(1), repmat([30; 100],1,sum(events_ID == events_unique(event_i))),"linestyle","--","marker","none","color",colors_event{event_i})
end
xlim([times(1) - days(5) times(end)])

% plot left hemisphere
nexttile(3,[2,1])
lfp_left_smov = filloutliers(lfp_left_s,"linear","movmedian",hours(24),"SamplePoints",time_left_s,"threshold",5);
timel =  time_left_s - time_left_s(1);
timel.Format = "d";
scatter(timel, lfp_left_smov,2,"k",'filled')

hold on
scatter(medchanges_date(medchanges < 0),4500 ,85,"kv","filled")
scatter(medchanges_date(medchanges > 0),4500,85,"k^","filled")
ylabel(" LFP Amplitude @12.7 Hz")
% ylim([0 2450])
ylim([0 5000])
xlim([timel(1) - days(5) timel(end)])
box off

% plot right hemisphere
nexttile(5,[2,1])
lfp_right_smov = filloutliers(lfp_right_s,"linear","movmedian",hours(24),"SamplePoints",time_right_s,"threshold",5);
timer =  time_right_s - time_right_s(1);
timer.Format = "d";
scatter(timer, lfp_right_smov,2,"k",'filled')
hold on
ylabel(" LFP Amplitude @7.81 Hz")
ylim([0 5000])
xlim([timer(1) - days(5) timer(end)])

if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_rawtrends_scattered"));
    saveas(gcf,figname,"png")
end

%% Plot clinical scales over time

figure("renderer","painters","Position",[400 400 750 300])
nexttile(1,[2,1])
scatter(sessions_clin,ybocs,40,"filled")
hold on
scatter(sessions_clin,madrs,40,"filled")
plot(sessions_clin,ybocs,"color",[0.8500, 0.3250, 0.0980],"linestyle","-")
plot(sessions_clin,madrs,"color",[0.8500, 0.3250, 0.0980],"linestyle","--")
box off
ylabel( " Clinical scale ")
%
if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_clintrends"));
    saveas(gcf,figname,"png")
end

%% Plot example of circadian rythm

circROI = [time_left_s(1) + days(154) time_left_s(1) + days(161)];

figure("renderer","painters","Position",[300   300   750 500])
t = tiledlayout(5,1);
% left hemisphere
nexttile(1,[2 1])
idx = isbetween(time_left_s,circROI(1),circROI(2));
timespan_circ_vec = time_left_s(idx);
hh = hms(timespan_circ_vec);
hh00 = find(hh == 0);
hh06 = find(hh == 6);
hh12 = find(hh == 12);
dayyy = time_left_s - time_left_s(1);
dayyy.Format = "d";
dayyy = dayyy(idx);
scatter(dayyy, zscore(lfp_left_smov(idx)),10,"k",'filled')
hold on
plot([dayyy(hh00); dayyy(hh00)], [-5 5],"color","k","linewidth",.2)
plot([dayyy(hh06); dayyy(hh06)], [-5 5],"color","r","linewidth",.2)
plot([dayyy(hh12); dayyy(hh12)], [-5 5],"color","b","linewidth",.2)
ylabel(" LFP Amplitude @12.7 Hz")
ylim([-5 5])
axis off
% right hemisphere
nexttile(3,[2 1])
scatter(dayyy, zscore(lfp_right_smov(idx)),10,"k",'filled')
hold on
plot([dayyy(hh00); dayyy(hh00)], [-5 5],"color","k","linewidth",.2)
plot([dayyy(hh06); dayyy(hh06)], [-5 5],"color","r","linewidth",.2)
plot([dayyy(hh12); dayyy(hh12)], [-5 5],"color","b","linewidth",.2)
ylabel(" LFP Amplitude @12.7 Hz")
ylim([-5 5])
axis off

if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_circ_example"));
    saveas(gcf,figname,"png")
end

%% Plot de-trended of circadian rythm
 
smooth = 24*6; % 24h moving median filter
 
lfp_left_smovd = movmedian(lfp_left_smov,smooth);
lfp_right_smovd = movmedian(lfp_right_smov,smooth); 
 
 
figure("renderer","painters","Position",[400   400   750 500])
t = tiledlayout(5,1);
nexttile(1,[1,1])
times =  time_left_s_stim - time_left_s(1);
times.Format = "d";
plot(times,stim_left_s,"r")
hold on
ylabel(" Stim Amp. [mA] ")
ylim([0 7])
xlim([times(1) - days(5) times(end)])
box off

nexttile(2,[2,1])
timel =  time_left_s - time_left_s(1);
timel.Format = "d";
scatter(timel, lfp_left_smovd,2,"k",'filled')
hold on
scatter(medchanges_date(medchanges < 0),4500 ,85,"kv","filled")
scatter(medchanges_date(medchanges > 0),4500,85,"k^","filled")
ylabel(" LFP Amplitude @12.7 Hz")
ylim([0 2450])
xlim([timer(1) - days(5) timer(end)])
box off
 
nexttile(4,[2,1])
timer =  time_right_s - time_right_s(1);
timer.Format = "d";
scatter(timer, lfp_right_smovd,2,"k",'filled')
hold on
ylabel(" LFP Amplitude @7.81 Hz")
ylim([0 2450])
xlim([timer(1) - days(5) timer(end)])
box off

if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_dailytrends"));
    saveas(gcf,figname,"png")
end

%% Focus initial LFP amplitude drop (S1 vs S3)

figure("renderer","painters","Position",[400   400   750 500])
t = tiledlayout(5,1);
nexttile(1,[2,1])
timel =  time_left_s - time_left_s(1);
timel.Format = "d";
scatter(timel, lfp_left_smovd,10,"k","filled")
hold on
xlim([timel(1) timel(3000)])
box off
 
nexttile(3,[2,1])
timer =  time_right_s - time_right_s(1);
timer.Format = "d";
scatter(timer, lfp_right_smovd,10,"k","filled")
hold on
xlim([timer(1) timer(3000)])
box off

if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_rawtrend_focuss"));
    saveas(gcf,figname,"png")
end

%% Circadian analysis
 
LFP = timetable(time_left_s',lfp_left_smov',lfp_right_smov');
LFP.Properties.VariableNames = {'E1L','E1R'};
newtimes = datetime('now') + hours(0:23);
 
Days = day(LFP.Time);
Hours = hour(LFP.Time);
Months = month(LFP.Time);
 
Days_end = [find(diff(Days) ~= 0); numel(Days)] ;
Days_start = [1; Days_end(1:end-1)+1];
Days_end= [Days_end(1:8); 1142; Days_end(9:end)];
Days_start= [Days_start(1:9); 1143; Days_start(10:end)];
 
Days_u = 1 : numel(Days_start);
Days = nan(numel(LFP.Time),1);
 
for day_i = 1 : numel(Days_start)
    Days(Days_start(day_i) : Days_end(day_i)) = day_i;
end
 
Days_unique = unique(Days);
Hours_unique = unique(Hours);
 
LFPcircday = nan(2,numel(Days_unique),numel(Hours_unique));
 
 
for hour_i = 1 : numel(Hours_unique)
    for day_i = 1 : numel(Days_unique)
        LFPcircday(1,day_i,hour_i) = nanmedian(LFP.E1R(Days == Days_unique(day_i) & Hours == Hours_unique(hour_i)));
        LFPcircday(2,day_i,hour_i) = nanmedian(LFP.E1L(Days == Days_unique(day_i) & Hours == Hours_unique(hour_i)));
    end
end
 
 
for day_i = 1: numel(Days_unique)
    tmp = squeeze(LFPcircday(1,day_i,:)/nansum(LFPcircday(1,day_i,:)));
    LFPcircH(1,day_i) = -sum(tmp.*log2(tmp));
    tmp = squeeze(LFPcircday(2,day_i,:)/nansum(LFPcircday(2,day_i,:)));
    LFPcircH(2,day_i) = -sum(tmp.*log2(tmp));
end
 
LFPcirc = zeros(2,numel(Hours_unique));
LFP25pth = zeros(2,numel(Hours_unique));
LFP75pth = zeros(2,numel(Hours_unique));
 
for hour_i = 1 : numel(Hours_unique)

    LFPcirc(1,hour_i) = nanmedian(LFP.E1R( Hours == Hours_unique(hour_i)));
    LFPcirc(2,hour_i) = nanmedian(LFP.E1L( Hours == Hours_unique(hour_i)));
    
    LFP25pth(1,hour_i) = squeeze(prctile(LFPcircday(1,:,hour_i),25));
    LFP25pth(2,hour_i) = squeeze(prctile(LFPcircday(2,:,hour_i),25));
    LFP75pth(1,hour_i) = squeeze(prctile(LFPcircday(1,:,hour_i),75));
    LFP75pth(2,hour_i) = squeeze(prctile(LFPcircday(2,:,hour_i),75));
end
LFPcirc_z = squeeze(nanmedian(zscore(LFPcircday,0,3),2));
LFP25pth_z = squeeze(prctile(zscore(LFPcircday,0,3),25,2));
LFP75pth_z = squeeze(prctile(zscore(LFPcircday,0,3),75,2));
 
 
figure("renderer","painters","position",[300 300 1200 400])
tiledlayout(1,3)
tmp = zscore(LFPcircday,0,3);
nexttile(2,[1,1])
hold on
for hour_i = 1:24
    scatter(hours(hour_i-1) + hours(0.1*randn(1,squeeze(numel(LFPcircday(1,:,hour_i))))),squeeze(tmp(1,:,hour_i)),'SizeData',1,"markerfacecolor",[102 204 0 ]/256,'markeredgecolor',[102 204 0 ]/256,'markerfacealpha',.01)
end
fill([0:hours(1):hours(23) fliplr(0:hours(1):hours(23))],[LFP75pth_z(1,:)  fliplr(LFP25pth_z(1,:))] ,[102 204 0 ]/256,"facealpha",.35,"linestyle","none")
plot(0:hours(1):hours(23),LFPcirc_z(1,:),"linewidth",1.5,"color",[102 204 0 ]/256)
xlabel(" Time [Hours] ")
ylabel(" LFP z-Amplitude @7.81 Hz")
box off
title("Right ")
ylim([ -3 3])
nexttile(1,[1,1])
hold on
for hour_i = 1:24
    scatter(hours(hour_i-1) + hours(0.1*randn(1,squeeze(numel(tmp(2,:,hour_i))))),squeeze(tmp(2,:,hour_i)),'SizeData',1,"markerfacecolor","b",'markeredgecolor',"b",'markerfacealpha',.01)
end
fill([0:hours(1):hours(23) fliplr(0:hours(1):hours(23))],[LFP75pth_z(2,:)  fliplr(LFP25pth_z(2,:))] ,"b","facealpha",.35,"linestyle","none")
plot(0:hours(1):hours(23),LFPcirc_z(2,:),"linewidth",1.5,"color","k")
xlabel(" Time [Hours] ")
ylabel(" LFP z-Amplitude @12.71 Hz")
title("Left ")
ylim([ -3 3])
box off
 
 
nexttile(3,[1,1])
hold on
tmp = zscore(LFPcircday,0,3);
 
for hour_i = 1:24
    scatter(hours(hour_i-1) + hours(0.1*randn(1,squeeze(numel(LFPcircday(2,:,hour_i))))),squeeze(mean(tmp(:,:,hour_i))),'SizeData',1,"markerfacecolor",[.4 .4 .4],'markeredgecolor',[.4 .4 .4],'markerfacealpha',.01)
end
fill([0:hours(1):hours(23) fliplr(0:hours(1):hours(23))],mean([[LFP75pth_z(1,:)  fliplr(LFP25pth_z(1,:))] ; [LFP75pth_z(2,:)  fliplr(LFP25pth_z(2,:))]]),[.8 .8 .8],"facealpha",.35,"linestyle","none")
plot(0:hours(1):hours(23),mean(LFPcirc_z),"linewidth",1.5,"color","k")
xlabel(" Time [Hours] ")
ylabel(" LFP z-Amplitude ")
title("Both ")
box off
ylim([-3 3])
 
if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_cyrcadian"));
    saveas(gcf,figname,"png")
end
 
%% Plot event-locked PSD markeb by the patient
 
 
figure("renderer","painters","position",[300 300 500 500])

sess_i = idx_forevent;
n_events = numel(events);
powEvents = struct();
 
 
Contact1_l = zeros(1,n_events);
Contact2_l = 2*ones(1,n_events);
Contact1_r = zeros(1,n_events);
Contact2_r = 2*ones(1,n_events);
 
idx_distcoord_l = 3*ones(1,n_events);
idx_distcoord_r = 3*ones(1,n_events);
 
 
% FOOOF settings
settings = struct();  % Use defaults
f_range = [1, 55];
settings.peak_width_limits = [0.5 12];
settings.max_n_peaks = 5;
settings.min_peak_height = .0;
settings.peak_threshold = 2;
settings.aperiodic_mode = 'fixed';
settings.verbose = true;
 
events_ID = [data_all{1, idx_forevent}.DiagnosticData.LfpFrequencySnapshotEvents.EventID];
events_names = {data_all{1, idx_forevent}.DiagnosticData.LfpFrequencySnapshotEvents.EventName};
for event_i = 1:n_events
    powEvents(event_i).Contact1_l = Contact1_l(event_i);
    powEvents(event_i).Contact2_l = Contact2_l(event_i);
    powEvents(event_i).Contact1_r = Contact1_l(event_i);
    powEvents(event_i).Contact2_r = Contact2_l(event_i);
    powEvents(event_i).freq = data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents(event_i).LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.Frequency;
    idx_Freq = dsearchn(powEvents(event_i).freq,[7.81 12.7]');
    powEvents(event_i).pow_left = data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents(event_i).LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.FFTBinData;
    powEvents(event_i).pow_right = data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents(event_i).LfpFrequencySnapshotEvents.HemisphereLocationDef_Right.FFTBinData;
    powEvents(event_i).date = data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents(event_i).DateTime;
    powEvents(event_i).eventID = data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents(event_i).EventID;
    powEvents(event_i).eventName = data_all{1, sess_i}.DiagnosticData.LfpFrequencySnapshotEvents(event_i).EventName;
    powEvents(event_i).pow_tracked_left = powEvents(event_i).pow_left(idx_Freq(2));
    powEvents(event_i).pow_tracked_right = powEvents(event_i).pow_right(idx_Freq(1));
    
    % Run FOOOF
    fooof_results = fooof(powEvents(event_i).freq', powEvents(event_i).pow_left', f_range, settings, true);
    fooof_plot(fooof_results)
    powEvents(event_i).pow_tracked_left_fooofed = fooof_results.power_spectrum(dsearchn(fooof_results.freqs',12.7)) - fooof_results.ap_fit(dsearchn(fooof_results.freqs',12.7))  ;
    powEvents(event_i).pow_left_fooofed = fooof_results.power_spectrum - fooof_results.ap_fit;
    powEvents(event_i).freq_fooof = fooof_results.freqs;
    fooof_results = fooof(powEvents(event_i).freq', powEvents(event_i).pow_right', f_range, settings, true);
    powEvents(event_i).pow_right_fooofed = fooof_results.power_spectrum - fooof_results.ap_fit;
    powEvents(event_i).pow_tracked_right_fooofed = fooof_results.power_spectrum(dsearchn(fooof_results.freqs',7.81)) - fooof_results.ap_fit(dsearchn(fooof_results.freqs',7.81))  ;
end
 
 
 
 
colors = {[0 0 0]};
figure("renderer","painters","Position",[400 400 800 700])
t = tiledlayout(2,2);
nexttile(1,[1,1])
hold on
for event_i = 1:numel(events_unique)
    
    idx_event = find(events_ID == events_unique(event_i));
    tmp = [];
    for plot_i = 1 : numel(idx_event)
        plot(powEvents(idx_event(plot_i)).freq,powEvents(idx_event(plot_i)).pow_left,"color",[colors{1,event_i} .4],"linewidth",.6)
        tmp = [tmp powEvents(idx_event(plot_i)).pow_left];
    end
    plot(powEvents(idx_event(plot_i)).freq,mean(tmp,2),"color",[colors{1,event_i}],"linewidth",1.5)
    plot([12.7 12.7],[0 5],"k--")
    
    
end
 
xlim([0 50])
xlabel(" Frequency [Hz] ")
ylabel(" PSD [\muV^{2}/Hz] ")
ylim([ 0 5])
 
nexttile(3,[1,1])
hold on
 
for event_i = 1:numel(events_unique)
    idx_event = find(events_ID == events_unique(event_i));
    tmp = [];
    for plot_i = 1 : numel(idx_event)
        plot(powEvents(idx_event(plot_i)).freq_fooof,10.^powEvents(idx_event(plot_i)).pow_left_fooofed,"color",[colors{1,event_i} .4],"linewidth",.6)
        tmp = [tmp;10.^powEvents(idx_event(plot_i)).pow_left_fooofed];
    end
    plot(powEvents(idx_event(plot_i)).freq_fooof,mean(tmp),"color",[colors{1,event_i}],"linewidth",1.5)
    plot([12.7 12.7],[0 3],"k--")
end
 
xlim([0 50])
xlabel(" Frequency [Hz] ")
ylabel(" PSD-fooofed [a.u.] ")
 
nexttile(2,[1,1])
hold on
 
for event_i = 1:numel(events_unique)
    idx_event = find(events_ID == events_unique(event_i));
    tmp = [];
    for plot_i = 1 : numel(idx_event)
        plot(powEvents(idx_event(plot_i)).freq,powEvents(idx_event(plot_i)).pow_right,"color",[colors{1,event_i} .4],"linewidth",.6)
        tmp = [tmp  powEvents(idx_event(plot_i)).pow_right];
    end
    plot(powEvents(idx_event(plot_i)).freq,mean(tmp,2),"color",[colors{1,event_i}],"linewidth",1.5)
    plot([7.81 7.81],[0 5],"k--")
end
 
xlim([0 50])
xlabel(" Frequency [Hz] ")
ylabel(" PSD [\muV^{2}/Hz] ")
ylim([ 0 5])
 
 
nexttile(4,[1,1])
hold on
 
for event_i = 1:numel(events_unique)
    idx_event = find(events_ID == events_unique(event_i));
    tmp = [];
    for plot_i = 1 : numel(idx_event)
        plot(powEvents(idx_event(plot_i)).freq_fooof,10.^powEvents(idx_event(plot_i)).pow_right_fooofed,"color",[colors{1,event_i} .4],"linewidth",.6)
        tmp = [tmp;10.^powEvents(idx_event(plot_i)).pow_right_fooofed];
    end
    plot(powEvents(idx_event(plot_i)).freq_fooof,mean(tmp),"color",[colors{1,event_i}],"linewidth",1.5)
    plot([7.81 7.81],[0 3],"k--")
end
 
xlim([0 50])
xlabel(" Frequency [Hz] ")
ylabel(" PSD-fooofed [a.u.] ")

if FLAG_FIGSAVE
    figname = fullfile(PATH_FIGURE,strcat(info_data.PATIENT,"_eventsPSD"));
    saveas(gcf,figname,"png")
end
 

 
 

