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


% -------- S2GC MODEL SET UP
% Defining the views (sets of features from the same matrix)
view1 = Data(:,'dispersion_nodal':'SHAPE_Volume')
ViewNum(1) = size(view1,2);
ViewIdx{1} = 1:ViewNum(1);

view2 = Data(:,'GLCM_Homogeneity':'GLCM_Dissimilarity')
ViewNum(2) = size(view2,2);
ViewIdx{2} = 1:ViewNum(2);

view3 = Data(:,'GLRLM_SRE':'GLRLM_RP')
ViewNum(3) = size(view3,2);
ViewIdx{3} = 1:ViewNum(3);

view4 = Data(:,'NGLDM_Coarseness':'NGLDM_Busyness')
ViewNum(4) = size(view4,2);
ViewIdx{4} = 1:ViewNum(4);

view5 = Data(:,'GLZLM_SZE':'GLZLM_ZP')
ViewNum(5) = size(view5,2);
ViewIdx{5} = 1:ViewNum(5);

% Process data as they are Cox-ready
X = Cox_data_processed(X(:,3:end),X(:,1),'censoring',~X(:,2)); % Whole dataset
X_hum = Cox_data_processed(X_hum(:,3:end),X_hum(:,1),'censoring', ...
    ~X_hum(:,2)); % ICH dataset
X_int = Cox_data_processed(X_int(:,3:end),X_int(:,1),'censoring', ...
    ~X_int(:,2)); % INT dataset

% Define S2GC model parameter as output by best_parameters.m
lambda = 0.5;
rho = 0.01;
tau = 0.1;
% OptimalClusterNum = 2;   

% % Displaying the final graph
% figure;
% Graph = graph(G);
% plot(Graph,'-.dr','NodeLabel',Linfoma_Patients)

% Perform S2GC model on the three datasets
[wm, funcValMv, G] = S2GC(X, lambda, rho, tau, ViewIdx); % Whole dataset
[wm_hum, funcValMv_hum, G_hum] = S2GC(X_hum, lambda, rho, ...
    tau, ViewIdx); % ICH dataset
[wm_int, funcValMv_int, G_int] = S2GC(X_int, lambda, rho, ...
    tau, ViewIdx); % INT dataset

% % Displaying the final graphs (whole, ICH, INT)
% figure;
% Graph = graph(G);
% plot(Graph,'-.dr','NodeLabel',Linfoma_Patients)

% figure;
% Graph_hum = graph(G_hum);
% plot(Graph_hum,'-.dr','NodeLabel',Patients_hum)

% figure;
% Graph_int = graph(G_int);
% plot(Graph_int,'-.dr','NodeLabel',Patients_int)

% Perform spectral clustering on the three datasets
OptimalClusterNum = find_k(G);
idx = spectralcluster(double(G),OptimalClusterNum, ...
    'Distance','precomputed','LaplacianNormalization','symmetric');
group = num2cell(num2str(idx));

OptimalClusterNum = find_k(G_hum);
idx_hum = spectralcluster(double(G_hum),OptimalClusterNum, ...
    'Distance','precomputed','LaplacianNormalization','symmetric');
group_hum = num2cell(num2str(idx_hum));

OptimalClusterNum = find_k(G_int);
% idx_int = SpectralClustering(double(G_int),OptimalClusterNum);
idx_int = spectralcluster(double(G_int),OptimalClusterNum, ...
    'Distance','precomputed','LaplacianNormalization','symmetric');
group_int = num2cell(num2str(idx_int));


% -------- RESULTS 
% Check survival curves differences of found groups
[p,fh,stats] = MatSurv(X.sorty,~X.cens,group);
[p_hum,fh_hum,stats_hum] = MatSurv(X_hum.sorty,~X_hum.cens,group_hum);
[p_int,fh_int,stats_int] = MatSurv(X_int.sorty,~X_int.cens,group_int);


% Estimate R2 of linear regression
mdl = fitlm(X.X,idx); % with radiomics variables
% fprintf('Ordinary R-squared for all dataset: %f \n',mdl.Rsquared.Ordinary);
fprintf('Adjusted R-squared for all dataset: %f \n', ...
    mdl.Rsquared.Adjusted);
fprintf('RMSE for all dataset: %f \n',mdl.RMSE);
fprintf('P-Value for all dataset: %d \n',anova(mdl,'summary').pValue(2));

mdl_hum = fitlm(X_hum.X,idx_hum);
% fprintf('Ordinary R-squared for HUM dataset: %f \n',mdl_hum.Rsquared.Ordinary);
fprintf('Adjusted R-squared for HUM dataset: %f \n', ...
    mdl_hum.Rsquared.Adjusted);
fprintf('RMSE for HUM dataset: %f \n',mdl_hum.RMSE);
fprintf('P-Value for HUM dataset: %d \n',anova(mdl_hum,'summary').pValue(2));

mdl_int = fitlm(X_int.X,idx_int);
% fprintf('Ordinary R-squared for INT dataset: %f \n',mdl_int.Rsquared.Ordinary);
fprintf('Adjusted R-squared for INT dataset: %f \n', ...
    mdl_int.Rsquared.Adjusted);
fprintf('RMSE for INT dataset: %f \n',mdl_int.RMSE);
fprintf('P-Value for INT dataset: %f \n',anova(mdl_int,'summary').pValue(2));

% Estimate pseudo R2 logistic regression
for index=1:length(idx)
    idx(idx==index) = index-1; % to make idx start from 0
end
mdl = fitglm(X.X,idx);
fprintf('Adjusted R-squared for all dataset: %f \n', ...
    mdl.Rsquared.Adjusted);

for index=1:length(idx_hum)
    idx_hum(idx_hum==index) = index-1;  % to make idx_hum start from 0
end
mdl_hum = fitglm(X_hum.X,idx_hum);
fprintf('Adjusted R-squared for HUM dataset: %f \n', ...
    mdl_hum.Rsquared.Adjusted);

for index=1:length(idx_int)
    idx_int(idx_int==index) = index-1;  % to make idx_int start from 0
end
mdl_int = fitglm(X_int.X,idx_int);
fprintf('Adjusted R-squared for INT dataset: %f \n', ...
    mdl_int.Rsquared.Adjusted);

% Perform two sided tests on groups
p_vect = [];
% On whole dataset
for iFeat=1:size(X.X,2)
    p_i = kruskalwallis(X.X(:,iFeat), group, 'off'); % more than 2 groups
    % p_i = ranksum(X.X(strcmp(group, '1'),iFeat), ...
    %     X.X(strcmp(group, '2'),iFeat), 'tail','both');  % only 2 groups
    p_vect = [p_vect p_i];
end

p_vect_hum = [];
% On ICH dataset
for iFeat=1:size(X_hum.X,2)
    p_i = kruskalwallis(X_hum.X(:,iFeat), group_hum, 'off'); % more than 2 groups
    % p_i = ranksum(X_hum.X(strcmp(group_hum, '1'),iFeat), ...
    %     X_hum.X(strcmp(group_hum, '2'),iFeat), 'tail','both');  % only 2 groups
    p_vect_hum = [p_vect_hum p_i];
end
p_vect_int = [];
% On INT dataset
for iFeat=1:size(X_int.X,2)
    p_i = kruskalwallis(X_int.X(:,iFeat), group_int, 'off'); % more than 2 groups
    % p_i = ranksum(X_int.X(strcmp(group_int, '1'),iFeat), ...
    %     X_int.X(strcmp(group_int, '2'),iFeat), 'tail','both');  % only 2 groups
    p_vect_int = [p_vect_int p_i];
end

threshold = 0.1; % 0.05
Linfoma_ColNames = Linfoma.Properties.VariableNames;


% -------- TRANSFER S2GC MODEL FROM ONE DOMAIN (ICH) TO THE OTHER (INT)
% Extract significant features for ICH
cols_boolean = p_vect_hum<threshold;
cols = Linfoma_ColNames(cols_boolean);

% Split training and testing 
X_train = X_hum.X(:,cols_boolean);
X_train_tbl = array2table(X_train, 'VariableNames',cols');
y_train = group_hum;
X_test = X_int.X(:,cols_boolean);
X_test_tbl = array2table(X_test, 'VariableNames',cols');
y_test = group_int;

% Fit the Tree ensemble
rf = TreeBagger(100, X_train_tbl, y_train, ...
    'Method', 'classification', 'OOBPrediction', 'on', ...
    'OOBPredictorImportance', 'on', 'NumPredictorsToSample', 'all', ...
    'MinLeafSize', 5, 'Prior', 'Empirical'); % 'Uniform'
% 'Cost' per ora è bilanciato
label = rf.predict(X_train);

disp('Traning performance:')
C_train = confusionmat(y_train,label)
acc_train = sum(diag(C_train))/sum(sum(C_train))
[statistics_train] = statsOfMeasure(C_train, 1);

% Check results
[p_train,fh_train,stats_train] = MatSurv(X_hum.sorty,~X_hum.cens,label);

y_pred = rf.predict(X_test_tbl);

disp('Testing performance:')
C_test = confusionmat(y_test,y_pred)
acc_test = sum(diag(C_test))/sum(sum(C_test))
[statistics_test] = statsOfMeasure(C_test, 1);

[p_test,fh_test,stats_test] = MatSurv(X_int.sorty,~X_int.cens,y_pred);

mdl_test = fitlm(X_int.X,str2double(y_pred));
% fprintf('Ordinary R-squared for INT dataset: %f \n',mdl_test.Rsquared.Ordinary);
fprintf('Adjusted R-squared for INT dataset: %f \n', ...
    mdl_test.Rsquared.Adjusted);
fprintf('RMSE for INT dataset: %f \n',mdl_test.RMSE);
fprintf('P-Value for INT dataset: %f \n',anova(mdl_test,'summary').pValue(2));

% Perform two sided tests on groups
p_vect_int_post = [];
for iFeat=1:size(X_int.X,2)
    % p_i = kruskalwallis(X_int.X(:,iFeat), y_pred, 'off');
    p_i = ranksum(X_int.X(strcmp(y_pred, '1'),iFeat), ...
        X_int.X(strcmp(y_pred, '2'),iFeat), 'tail','both');
    p_vect_int_post = [p_vect_int_post p_i];
end

cols_boolean_post = p_vect_int_post<threshold;
cols_post = Linfoma_ColNames(cols_boolean_post);


% -------- EXPLAINABILITY
% Importance plot
imp_ = rf.OOBPermutedPredictorDeltaError;
[imp,sort_indexing] = sort(imp_, 'descend');
imp = imp(abs(imp)>0.2);
importance_names_ = cols(sort_indexing); % rf.PredictorNames(sort_indexing);
importance_names = importance_names_(abs(imp)>0.15);
for iVar=1:length(importance_names)
    var = importance_names(iVar);
    importance_names(iVar) = strrep(var,'_','\_');
end
importance_names = categorical(importance_names);
importance_names = reordercats(importance_names,string(importance_names));

figure;
barh(importance_names, imp, 'FaceColor', '#A9A9A9'); % black
title('Rule Extraction Model - Importance Plot');
xlabel('Predictor importance estimates');
ylabel('Predictors');

% Rule extraction: see R script



%%
function k = find_k(G)

    OptimalClusterNum = 10; % big enough

    [idx,V,D] = spectralcluster(double(G),OptimalClusterNum, ...
        'Distance','precomputed','LaplacianNormalization','randomwalk');

    zeros_eig = find(D==0);
    k = zeros_eig(end)+1;

    figure;
    scatter(1:OptimalClusterNum,D)
 
end
