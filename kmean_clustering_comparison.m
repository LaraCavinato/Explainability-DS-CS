close all
clear
clc

% Add folders and subfolders to the directory
addpath(genpath(pwd))

% -------- DATASETS PREPROCESSING
% Load the datasets in DATA/ (not included in this repo)
Linfoma_hum_ = readtable("DATA/ich_linfoma.xls");
Linfoma_int = readtable("DATA/int_linfoma.xls");

Linfoma = [Linfoma_hum; Linfoma_int]; % merging the two

num_hum_patients = size(Linfoma_hum,1); % # of patients in ICH dataset
num_int_patients = size(Linfoma_int,1); % # of patients in INT dataset

Linfoma_Patients = Linfoma.Pt_ID; % Patients' IDs
Response = Linfoma(:,end-1:end); % Responde variable
Linfoma = Linfoma(:,2:end-4); % Radiomic features

Linfoma.AGE = int(round(Linfoma.AGE)); % Casting AGE to integer

Data = table2array(Linfoma); % Table into matrix
Data = Standard_Normalization(Data); % Normalize data
Data = standardizeCols(Data); % Standardize data

% Dividing the two datasets
Data_hum = Data(1:num_hum_patients,:);
Data_int = Data(num_hum_patients+1:num_hum_patients+num_int_patients,:);

% Getting survival response for the two datasets
Response_hum = Response(1:num_hum_patients,:);
Response_int = Response(num_hum_patients+1:num_hum_patients+ ...
    num_int_patients,:);

% Getting patients IDs for the two datasets
Patients_hum = Linfoma_Patients(1:num_hum_patients);
Patients_int = Linfoma_Patients(num_hum_patients+1:num_hum_patients+ ...
    num_int_patients);

% Final data to be used
y = table2array(Response);
X = [y Data];

y_hum = table2array(Response_hum);
y_int = table2array(Response_int);
X_hum = [y_hum Data_hum];
X_int = [y_int Data_int];


% -------- K-MEANS CLUSTERING SET UP
idx_kmeans = kmeans(X.X,2,'EmptyAction','drop');
group_kmeans = num2cell(num2str(idx_kmeans));
% Display results
[p_kmeans,fh_kmeans,stats_kmeans] = MatSurv(X.sorty,~X.cens, ...
    group_kmeans);

idx_kmeans_hum = kmeans(X_hum.X,2,'EmptyAction','drop');
group_kmeans_hum = num2cell(num2str(idx_kmeans_hum));
% Display results
[p_kmeans_hum,fh_kmeans_hum,stats_kmeans_hum] = MatSurv(X_hum.sorty ...
    ,~X_hum.cens,group_kmeans_hum);

idx_kmeans_int = kmeans(X_int.X,2,'EmptyAction','drop');
group_kmeans_int = num2cell(num2str(idx_kmeans_int));
% Display results
[pkmeans_int,fhkmeans_int,statskmeans_int] = MatSurv(X_int.sorty, ...
    ~X_int.cens,group_kmeans_int);
