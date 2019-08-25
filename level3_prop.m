function [m_p,pr_p]=level3_prop(dpm)
global C prop_pr_forest prop_pr_DOE MidPara DsPm

%% mass estimation

    dpi=round(dpm/0.0254);      % convert diameter unit from metre to inch
    m_p=(0.0036907*dpi^3-0.021325*dpi^2+0.9344*dpi)/1000;
    
    % for large propeller, the price shows a linear relationship, while for
    % smaller propellers, it's more like a poly3 relationship. the two
    % correlation cross at dpi=13.11, so 13 inch is the turning point
%% price estimation
% fprintf('Predicting Validation set\n');
No_Trees_pr=length(prop_pr_forest.Trees);
X=dpi;
for i=1:No_Trees_pr
%     fprintf('Tree #%d\n',i);
    % define the leaf ID
    [~,leaf_ID] = predict(prop_pr_forest.Trees{i},X);
    
    for j=1:prop_pr_DOE{i}.NO_LEAF
        if sum(prop_pr_DOE{i}.LEAF_IDS(j)==leaf_ID)>0
%             fprintf('\tLeaf #%d\n',j);
            % calcualte the prediction
            [Y_PRED,SD_PRED] = predict(prop_pr_DOE{i}.LEAF{j}.MODEL,X(leaf_ID==prop_pr_DOE{i}.LEAF_IDS(j),:));            
            % store the predictions
            Y_p(leaf_ID==prop_pr_DOE{i}.LEAF_IDS(j),i) = Y_PRED;
            SD_p(leaf_ID==prop_pr_DOE{i}.LEAF_IDS(j),i) = SD_PRED;
           
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
pr_p=Price_p;  %!!!!!needs to be finished

end