% Code for running T. Helble TDOA tool on DCLDE 2013 train dataset

close all; clear all; clc;
cd('/home/kpalmer/AnacondaProjects/Localisation/Scripts');
%load('DCLDE2013_RW_localizations_NOPP6_EST__all14_timed1_95.mat')
load('Cornell_Cape_Cod_Setup.mat')
addpath('/home/kpalmer/AnacondaProjects/Localisation/Scripts');
addpath('/home/kpalmer/AnacondaProjects/Localisation/Localization_v1')
% Create the parameter file
%parm_dclde2013 = GPL_parameter_input


parm.plot =0; % Turn to 1 to see graphical represenation of detections
parm.pad =0; %defines overlap window in seconds for GPL to detect across boundaries, set to 0 for this project.
parm.loop =1; %Defaults to 5 druing setup process, should remain 1



% Load DCLDE meta data
% Load/Create Hydrophone Structure
dclde_2013_meta = xlsread('/cache/kpalmer/quick_ssd/data/DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');

% Create labels for plotting locations (undocumented function)
labels=sprintfc('%d',1:10);

plot(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 'o')
text(dclde_2013_meta(:,12), dclde_2013_meta(:,11),labels,'VerticalAlignment','top','HorizontalAlignment','left')

% Convert the meta data to a structure in the format that the
% GPL/localisation code expects

hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
    
end


% Create the array structure
array_struct.master=5; %choose your master hydrophone (in relation to definition in hydrophone_struct)
array_struct.slave=[1 2 3 6 7 8 9 10]; %define supporting hydrophones (ch 4 left out b/c it didn't produce data)
array_struct.latgrid=[41.904:.002: 42.42];  %defines total search area for localizations, go beyond array, can always truncate later
array_struct.longrid=[-70.58:.002:-70.027];
plot(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 'o')
text(dclde_2013_meta(:,12), dclde_2013_meta(:,11),labels,'VerticalAlignment','top','HorizontalAlignment','left')
xlim(minmax(array_struct.longrid))
ylim(minmax(array_struct.latgrid))


[array_struct] = setup_TDOA_grid(hydrophone_struct, array_struct, localize_struct, parm);

%[localize_struct,array_struct]=localize_parameter_input(hydrophone_struct,array_struct,parm);




localize_struct.parm.ssp_choice = 1; %use constant sound speed
localize_struct.parm.grid_depth = 5; %approx expected depth of whale
localize_struct.parm.num_calls = 1; %cross correlate one call at a time (standard)
localize_struct.parm.number_maxima = 1; %take only the maximum cross_corr peak
localize_struct.parm.lsq_cutoff = .1 %unused, later can filter on LSQ
localize_struct.parm.pow = 1.5 %FFT power of the call to cross corelate
%localize_struct.parm.ssp = ([depth (pos), ssp]); not used here but defines ssp for deeper water use.
%localize_struct.parm.phase=600;




localize_struct.parm.min_pairs = 2; % Set how many hyrophone pairs in solution (for a 6 hydrophone deployment, usually use value of 4 or 5).
% to save time, do this once and save. computes the expected time delays on all channels.




%% Run GPL

num_slates_file=15;
%[filename, pathname]= uigetfile('*.aif');
cwd='/cache/kpalmer/quick_ssd/data/DCLDE_2013_10Channel/All';
cd(cwd)
file_list=dir('*.aif');
cd(cwd);

prefix = 'DCLDE_2013_10_Chan';

workpath=pwd;
%addpath(pathname);
%cd(pathname);
cd(workpath);

gstart=1;

finish=length(file_list)-1;


for(j=1:length(hydrophone_struct))
    hyd(j).detection.calls=[];
    hyd(j).detection.parm=parm;
end
k=0;


extract_date=file_list(gstart).name(11:end-4);
start_date= datenum(str2num(extract_date(1:4)),...
    str2num(extract_date(5:6)),...
    str2num(extract_date(7:8)),...
    str2num(extract_date(10:11)),...
    str2num(extract_date(12:13)),...
    str2num(extract_date(13:14)));

%%
%for file_index=gstart:2

for file_index=gstart:finish
    
    % Report percentage of the job done on the screen
    perccount(file_index,finish);
    
    k=k+1;
    tic
    file_list(file_index).name;

    % Load the aiff file
    data=audioread(file_list(file_index).name)';

    
   
    % Iterate through the number of hydrophones
    for j=1:length(hydrophone_struct)
        
        if j ~=4
            [calls]=process_RW_sub(data(hydrophone_struct(j).channel,:),...
                parm,...
                start_date+datenum(0,0,0,0,0,parm.nrec/parm.sample_freq*num_slates_file*(k-1)),...
                num_slates_file);
            
            [sz1 sz2]=size(calls);
            
            for(jj=1:sz2)
                calls(jj).start_time=calls(jj).start_time+((k-1)*parm.nrec*num_slates_file);
                calls(jj).end_time=calls(jj).end_time+((k-1)*parm.nrec*num_slates_file);
            end
            hyd(j).detection.calls=[hyd(j).detection.calls,calls];
        end
        
    end
   
 toc
    
end


sizeOfSave=whos('calls');

if sizeOfSave.bytes <= 2000000000
    save(strcat('DCLDE2013_RW_localizations_',prefix,'_all14_timed_oct31_calls',num2str(gstart),'_',num2str(finish)),'calls')
else save(strcat('DCLDE2013_RW_localizations_',prefix,'_all14_timed_oct31_calls',num2str(gstart),'_',num2str(finish)),'calls', '-v7.3')
end

%% 

% Load the Calls 


% Run 2d localsation in loop
all_hyd = [1 2 3 5 6 7 8 9 10];



% Loop through all the hydrophones and create localisations/tdoa spaces

%for ii=1:length(all_hyd)
for ii=1:length(all_hyd)
    
    % Create the array structure
    array_struct.master= all_hyd(ii); 
    
    % Knock out the master hydrophone from the list of slaves
    new_hyd_array = all_hyd;
    new_hyd_array(ii) =[];
    
    array_struct.slave=new_hyd_array; %define supporting hydrophones (ch 4 left out b/c it didn't produce data)

    [array_struct] = setup_TDOA_grid(hydrophone_struct, array_struct,...
        localize_struct, parm);
    [localize_struct] = localize_cross_corr_index_check(array_struct,hyd,...
        localize_struct,all_hyd(ii));
    [localize_struct] = localize_LSQ_2D(array_struct, hydrophone_struct,hyd,...
        localize_struct,all_hyd(ii));
   
    % Add the hydrophone specific localization structure
    localize_struct.hyd(all_hyd(ii)).array_struct = array_struct;

    
end
% Save everything 

sizeOfSave=whos('localize_struct');
if sizeOfSave.bytes <= 2000000000
    
    save(strcat('DCLDE2013_RW_localizations',prefix,'_all_',num2str(gstart),'_',num2str(finish)),'localize_struct')
else save(strcat('DCLDE2013RW_localizations',prefix,'_all_',num2str(gstart),'_',num2str(finish)),'localize_struct', '-v7.3')
end



clear data            %swm added Sep 25, 2014 to ensure not impacting size of save
sizeOfSave=whos('hyd');

if sizeOfSave.bytes <= 2000000000
    save(strcat('DCLDE2013_RW_localizations_',prefix,'_all14_timed_oct31_hyd',num2str(gstart),'_',num2str(finish)),'hyd')
else save(strcat('DCLDE2013_RW_localizations_',prefix,'_all14_timed_oct31_hyd',num2str(gstart),'_',num2str(finish)),'hyd', '-v7.3')
end







%% 2 localisation (??)
sizeOfSave=whos('localize_struct');
if sizeOfSave.bytes <= 2000000000
    save(strcat('DCLDE2013_RW_localizations_',prefix,'_all14_timed',num2str(gstart),'_',num2str(finish)),'localize_struct')
else save(strcat('DCLDE2013RW_localizations_',prefix,'_all14_timed',num2str(gstart),'_',num2str(finish)),'localize_struct', '-v7.3')
end


%NOTE USING 2D solution
[localize_struct] = localize_LSQ_2D(array_struct, hydrophone_struct,hyd,localize_struct,1);



sizeOfSave=whos('localize_struct');
if sizeOfSave.bytes <= 2000000000
    save(strcat('RW_localizations_',prefix,'_all14_timed',num2str(gstart),'_',num2str(finish)),'localize_struct')
else save(strcat('RW_localizations_',prefix,'_all14_timed',num2str(gstart),'_',num2str(finish)),'localize_struct', '-v7.3')
end


%% Plot Results

min_lsq_score=.01;
sample_freq=2000;
start_time=0;


for(jj=1:5)
    
    figure(jj);
    for(j=1:length(hydrophone_struct))
        plot(hydrophone_struct(j).location(1,2),hydrophone_struct(j).location(1,1),'k*');
        hold on;
    end
    
    
    
    for j=1:length(localize_struct.hyd)
        hyd_number=j;
        
        if (isempty(localize_struct.hyd(hyd_number).coordinates)==0)
            
            score=localize_struct.hyd(hyd_number).score(jj,:);
            [k1, k2]=find(score < min_lsq_score);
            
            coordinates=localize_struct.hyd(hyd_number).coordinates;
            
            test=ndims(coordinates);
            
            if(test==3)
                coordinates=coordinates(jj,:,:);
                coordinates=squeeze(coordinates);
            end
            
            times=localize_struct.hyd(hyd_number).rtimes;
            
            
             figure(jj);
            hold on
            scatter(coordinates(2,k2), coordinates(1,k2), 20, times(k2)/sample_freq/3600, 'filled');
            
            
            time1900 = datenum(0,0,0,0,0,times(k2)/sample_freq)+start_time - datenum('1900-01-01','yyyy-mm-dd');
            time1900 = time1900+2;
            
            %  combined=[time1900; coordinates(1,k2);coordinates(2,k2); score(k2); ones(size(k2))*5; ones(size(k2))*1; ones(size(k2))*1; ones(size(k2))*1; ones(size(k2))*1];
            %  dlmwrite(strcat(name,'_viewer.dat'),combined','newline','pc','precision',15, '-append','delimiter', '\t');
            
        end
    end
end





























