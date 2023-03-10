function [n_split, n_trees, lr] = findBestEnsTree(X,Y)

% Cross-validate a deep classification tree and a stump. 
% These classification trees serve as benchmarks.

MdlDeep = fitctree(X,Y,'CrossVal','on','MergeLeaves','off', ...
    'MinParentSize',1);
MdlStump = fitctree(X,Y,'MaxNumSplits',1,'CrossVal','on');

% Cross-validate an ensemble of 150 boosted classification 
% trees using 5-fold cross-validation.
% Maximum number of splits no grater than n-1

n = size(X,1);
m = floor(log(n - 1)/log(3));
learnRate = [0.1 0.25 0.5 1];
numLR = numel(learnRate);
maxNumSplits = 3.^(0:m);
numMNS = numel(maxNumSplits);
numTrees = 150;
Mdl = cell(numMNS,numLR);

for k = 1:numLR
    for j = 1:numMNS
        t = templateTree('MaxNumSplits',maxNumSplits(j));
        Mdl{j,k} = fitcensemble(X,Y,'NumLearningCycles',numTrees,...
            'Learners',t,'KFold',5,'LearnRate',learnRate(k));
    end
end

% Estimate cumulative cross-validated misclassification error

kflAll = @(x)kfoldLoss(x,'Mode','cumulative');
errorCell = cellfun(kflAll,Mdl,'Uniform',false);
error = reshape(cell2mat(errorCell),[numTrees numel(maxNumSplits) ...
    numel(learnRate)]);
errorDeep = kfoldLoss(MdlDeep);
errorStump = kfoldLoss(MdlStump);

% Plot performance

mnsPlot = [1 round(numel(maxNumSplits)/2) numel(maxNumSplits)];
figure
for k = 1:3
    subplot(2,2,k)
    plot(squeeze(error(:,mnsPlot(k),:)),'LineWidth',2)
    axis tight
    hold on
    h = gca;
    plot(h.XLim,[errorDeep errorDeep],'-.b','LineWidth',2)
    plot(h.XLim,[errorStump errorStump],'-.r','LineWidth',2)
    plot(h.XLim,min(min(error(:,mnsPlot(k),:))).*[1 1],'--k')
    h.YLim = [0 0.2];    
    xlabel('Number of trees')
    ylabel('Cross-validated misclass. rate')
    title(sprintf('MaxNumSplits = %0.3g', maxNumSplits(mnsPlot(k))))
    hold off
end
hL = legend([cellstr(num2str(learnRate','Learning Rate = %0.2f')); ...
        'Deep Tree';'Stump';'Min. misclass. rate']);
hL.Position(1) = 0.6;

% Identify the maximum number of splits, number of trees, 
% and learning rate that yields the lowest misclassification 
% rate overall.

[minErr,minErrIdxLin] = min(error(:));
[idxNumTrees,idxMNS,idxLR] = ind2sub(size(error),minErrIdxLin);

% fprintf('\nMin. misclass. rate = %0.5f',minErr)
% fprintf('\nOptimal Parameter Values:\nNum. Trees = %d',idxNumTrees);
% fprintf('\nMaxNumSplits = %d\nLearning Rate = %0.2f\n',...
%     maxNumSplits(idxMNS),learnRate(idxLR))

n_split = maxNumSplits(idxMNS);
n_trees = idxNumTrees;
lr = learnRate(idxLR);

end
