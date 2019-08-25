function [m_motor,pr_motor]=level3_motor(param_motor)

global C Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest MidPara DsPm
% read the forest models, the forest including 13 trees, each leave on the
% tree is a Gaussian Regresison Process model. 
% Given the motor parameters, this function works out an estimated motor 
% mass for each tree, then select the one that has least error to the
% average motor/MTOM ratio. 
    
X=param_motor(1:3);     % read motor performance parameters
No_Trees_m=length(Motor_forest.Trees);
No_Trees_pr=length(Mdl_pr_forest.Trees);
%Y_p = nan(length(MASS_T),No_Trees);


%% mass estimation
%fprintf('Predicting Validation set\n');
for i=1:No_Trees_m
    %fprintf('Tree #%d\n',i);
    % define the leaf ID
    [~,leaf_ID] = predict(Motor_forest.Trees{i},X);
    
    for j=1:DOE_forest{i}.NO_LEAF
        if sum(DOE_forest{i}.LEAF_IDS(j)==leaf_ID)>0
            %fprintf('\tLeaf #%d\n',j);
            % calcualte the prediction
            [Y_PRED,~] = predict(DOE_forest{i}.LEAF{j}.MODEL,X(leaf_ID==DOE_forest{i}.LEAF_IDS(j),:));            
            % store the predictions
            Y_p(leaf_ID==DOE_forest{i}.LEAF_IDS(j),i) = Y_PRED/1000;
            %SD_p(leaf_ID==DOE_forest{i}.LEAF_IDS(j),i) = SD_PRED;
           
        end
    end
   % clc
end

%fprintf('Predicting based on minimum ratio error method:\n');
rMTOM=Y_p.*DsPm(1)./MidPara(5,1);      % work out the motor/MTOM 
% ratio for each estimated motor mass. Here takes the MTOM rather than
% GTOM because the correlation is between the motor mass and the MTOM,
% it reflects a rough ability of lifting load.
eMTOM=abs(rMTOM-C(30));     % work out the ratio error with an average value
[~, p_e]=min(eMTOM');       % locate the minimum error tree

m_motor=Y_p(p_e);           % select the mass estimation from this tree.


%% price estimation
% fprintf('Predicting Validation set\n');
for i=1:No_Trees_pr
%     fprintf('Tree #%d\n',i);
    % define the leaf ID
    [~,leaf_ID] = predict(Mdl_pr_forest.Trees{i},X);
    
    for j=1:DOE_pr_forest{i}.NO_LEAF
        if sum(DOE_pr_forest{i}.LEAF_IDS(j)==leaf_ID)>0
%             fprintf('\tLeaf #%d\n',j);
            % calcualte the prediction
            [Y_PRED,SD_PRED] = predict(DOE_pr_forest{i}.LEAF{j}.MODEL,X(leaf_ID==DOE_pr_forest{i}.LEAF_IDS(j),:));            
            % store the predictions
            Y_p(leaf_ID==DOE_pr_forest{i}.LEAF_IDS(j),i) = Y_PRED;
            SD_p(leaf_ID==DOE_pr_forest{i}.LEAF_IDS(j),i) = SD_PRED;
           
        end
    end
   % clc
end

% fprintf('Predicting based on weighted combination method:\n');
for i=1:size(Y_p,1)   
    % define the matrix C0
    C0 = inv(diag(SD_p(i,:)));
    W (i,:)= (C0*ones(No_Trees_pr,1))./(ones(1,No_Trees_pr)*C0*ones(No_Trees_pr,1));
    Price_p(i,1) = sum(Y_p(i,:).*W(i,:));    
end
pr_motor=Price_p;  %!!!!!needs to be finished
end