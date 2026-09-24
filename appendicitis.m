clc;
clear;
close all;

%% 1. 读取7个专家输入数据

file_path = 'D:\Desktop\阑尾炎数据\30shiyan.xlsx';
% file_path = 'D:\Desktop\阑尾炎数据\appendicitis_7groups.xlsx';



duixiang = 106;
shuxing = 7;
k_group = 7;

%% 2. 真实决策类别

POS0 = [zeros(1,21), ones(1,85)];
NEG0 = [ones(1,21), zeros(1,85)];

%% 3. 构造7个专家的数据

input_data = cell(1,k_group);

for p = 1:k_group
    sheet_name = sprintf('Sheet%d',p);
    Xp = readmatrix(file_path,'Sheet',sheet_name,'Range','A1:G106');
    Xtemp = nan(duixiang,shuxing);
    nr = min(size(Xp,1),duixiang);
    nc = min(size(Xp,2),shuxing);

    if nr > 0 && nc > 0
        Xtemp(1:nr,1:nc) = Xp(1:nr,1:nc);
    end

    input_data{p} = Xtemp;
end

%% 4. 相对损失参数

rho = 0.25;

%% 5. 建立各专家有效对象掩码

valid_mask = cell(1,k_group);

for p = 1:k_group
    Xp = input_data{p};
    valid_mask{p} = all(isfinite(Xp),2);
end

%% 6. 初始化变量

POS = zeros(k_group,duixiang);
BND = zeros(k_group,duixiang);
NEG = zeros(k_group,duixiang);

W_num = nan(k_group,shuxing);
PP = cell(1,k_group);
IDX_ALL = cell(1,k_group);
S_ALL = cell(1,k_group);
D_STATE = false(k_group,duixiang);
C_DOM_ALL = cell(1,k_group);
ALPHA_ALL = cell(1,k_group);
BETA_ALL = cell(1,k_group);

%% 7. 对每一个专家分别进行三支决策分类

for p = 1:k_group

    Xfull = input_data{p};
    validRowMask = valid_mask{p};
    validIdx = find(validRowMask);
    nValid = numel(validIdx);

    if nValid < 2
        error('专家%d有效对象数量不足2，无法继续计算。',p);
    end

    Pro_full = nan(duixiang,1);
    idx_full = nan(duixiang,1);
    xbar_full = nan(duixiang,1);
    alpha_full = nan(duixiang,1);
    beta_full = nan(duixiang,1);

    X = Xfull(validRowMask,1:shuxing);

    if any(~isfinite(X(:)))
        error('专家%d进行马氏距离K-means聚类前，X中包含NaN或Inf。',p);
    end

    %% 7.1 熵权法 Eq.(8) 到 Eq.(10)

    sumX = sum(X,1);

    if any(sumX <= 0)
        badCols = find(sumX <= 0);
        error('专家%d的Eq.(8)无法计算，属性列%s的列和小于等于0。',p,mat2str(badCols));
    end

    P = X ./ sumX;

    plnp = zeros(size(P));
    positiveMaskP = P > 0;
    plnp(positiveMaskP) = P(positiveMaskP) .* log(P(positiveMaskP));

    entropy_k = 1 / log(nValid);
    H = -entropy_k * sum(plnp,1);

    d = 1 - H;
    weightDen = sum(d);

    if weightDen == 0
        error('专家%d的Eq.(10)权重分母为0，论文未定义该退化情形。',p);
    end

    W = d / weightDen;
    W_num(p,:) = W;

    %% 7.2 加权相对优势关系 Eq.(11) 到 Eq.(15)

    C_dom = zeros(nValid,nValid);

    for i = 1:nValid
        for r = 1:nValid
            diff_ir = X(i,:) - X(r,:);
            A_ir = sum(W .* max(diff_ir,0));
            L_ir = sum(W .* max(-diff_ir,0));
            H_ir = A_ir + L_ir;

            if H_ir > eps
                C_dom(i,r) = A_ir / H_ir;
            else
                C_dom(i,r) = 0.5;
            end
        end
    end

    C_DOM_ALL{p} = C_dom;

    %% 7.3 马氏距离K-means聚类 Eq.(16)

    k_class = 2;

    options = struct();
    options.useBiasCov = false;
    options.MaxIter = 200;
    options.Replicates = 100;
    options.Start = 'plus';
    options.Seed = 1;
    options.reg = 1e-6;

    [idx,L] = kmeans_mahal_global(X,k_class,options);

    %% 7.4 状态集合 Eq.(17)

    center_score = L * W';
    [~,state_class] = max(center_score);
    D_valid = idx == state_class;
    D_STATE(p,validIdx(D_valid)) = true;

    %% 7.5 条件概率 Eq.(18)

    Pro_valid = zeros(nValid,1);

    for i = 1:nValid
        dominance_class = find(C_dom(:,i) >= 0.5);

        if isempty(dominance_class)
            Pro_valid(i) = 0;
        else
            Pro_valid(i) = sum(D_valid(dominance_class)) / numel(dominance_class);
        end
    end

    Pro_full(validRowMask) = Pro_valid;
    idx_full(validRowMask) = idx;

    PP{p} = Pro_full';
    IDX_ALL{p} = idx_full;

    %% 7.6 相对损失函数和决策阈值 Eq.(19)

    xbar = X * W';
    xbar_max = max(xbar);
    xbar_min = min(xbar);

    alpha = zeros(nValid,1);
    beta = zeros(nValid,1);

    for i = 1:nValid
        A = xbar_max - xbar(i);
        B = xbar(i) - xbar_min;

        den_alpha = (1-rho) * A + rho * B;
        den_beta = rho * A + (1-rho) * B;

        if den_alpha > eps
            alpha(i) = ((1-rho) * A) / den_alpha;
        else
            alpha(i) = 0;
        end

        if den_beta > eps
            beta(i) = (rho * A) / den_beta;
        else
            beta(i) = 0;
        end
    end

    xbar_full(validRowMask) = xbar;
    alpha_full(validRowMask) = alpha;
    beta_full(validRowMask) = beta;

    S_ALL{p} = xbar_full;
    ALPHA_ALL{p} = alpha_full;
    BETA_ALL{p} = beta_full;

    %% 7.7 三支决策 Eq.(20)

    for ii = 1:nValid
        obj = validIdx(ii);
        prob = Pro_valid(ii);

        if prob >= alpha(ii)
            POS(p,obj) = 1;
        elseif prob > beta(ii) && prob < alpha(ii)
            BND(p,obj) = 1;
        else
            NEG(p,obj) = 1;
        end
    end
end

%% 8. 专家信任系数 Eq.(22)

theta_expert = zeros(1,k_group);

for p = 1:k_group
    theta_expert(p) = sum(valid_mask{p}) / duixiang;
end

%% 9. 初始共识集合 Eq.(23)

weighted_POS0 = theta_expert * POS;
weighted_BND0 = theta_expert * BND;
weighted_NEG0 = theta_expert * NEG;

consensus_tol = 1e-12;

POS_T = zeros(1,duixiang);
BND_T = zeros(1,duixiang);
NEG_T = zeros(1,duixiang);

max_POS0 = max(weighted_POS0);

if max_POS0 > consensus_tol
    POS_T(abs(weighted_POS0 - max_POS0) < consensus_tol) = 1;
end

max_BND0 = max(weighted_BND0);

if max_BND0 > consensus_tol
    BND_T(abs(weighted_BND0 - max_BND0) < consensus_tol) = 1;
end

max_NEG0 = max(weighted_NEG0);

if max_NEG0 > consensus_tol
    NEG_T(abs(weighted_NEG0 - max_NEG0) < consensus_tol) = 1;
end

POS_consensus = POS_T;
BND_consensus = BND_T;
NEG_consensus = NEG_T;

%% 10. 迭代冲突消解 Eq.(25) 到 Eq.(36)

Su = zeros(1,k_group);
Influence = zeros(1,k_group);
Do = zeros(1,k_group);

Co = nan(duixiang,1);

iteration = 0;
resolved = false(duixiang,1);
max_iterations = duixiang + 1;

while any(~resolved) && iteration < max_iterations

    iteration = iteration + 1;
    unresolved_indices = find(~resolved);

    if isempty(unresolved_indices)
        break;
    end

    Su = calc_support(POS,BND,NEG,POS_T,BND_T,NEG_T,theta_expert,k_group);

    Influence = calc_influence_unresolved(POS,BND,NEG,theta_expert,resolved,valid_mask,k_group);

    Do = Su .* Influence;

    Co = nan(duixiang,1);

    for ii = 1:numel(unresolved_indices)

        i = unresolved_indices(ii);

        weighted_P = sum(Do .* POS(:,i)');
        weighted_B = sum(Do .* BND(:,i)');
        weighted_N = sum(Do .* NEG(:,i)');

        W_region = [weighted_P,weighted_B,weighted_N];
        total_weight = sum(W_region);

        if total_weight <= eps
            phi = [1/3,1/3,1/3];
        else
            phi = W_region / total_weight;
        end

        phi_nonzero = phi(phi > 0);

        Co(i) = -sum(phi_nonzero .* log(phi_nonzero)) / log(3);
    end

    Co_unresolved = Co(unresolved_indices);

    if any(~isfinite(Co_unresolved))
        error('第%d次迭代过程中出现非有限冲突度。',iteration);
    end

    min_co_val = min(Co_unresolved);
    co_tolerance = 1e-10;

    target_objs = unresolved_indices(abs(Co_unresolved - min_co_val) < co_tolerance);

    if isempty(target_objs)
        error('冲突消解过程中未能选择待处理对象。');
    end

    for t = 1:numel(target_objs)

        obj_idx = target_objs(t);

        vote_POS = sum(Do .* POS(:,obj_idx)');
        vote_BND = sum(Do .* BND(:,obj_idx)');
        vote_NEG = sum(Do .* NEG(:,obj_idx)');

        votes = [vote_POS,vote_BND,vote_NEG];
        max_vote = max(votes);
        vote_tolerance = 1e-10;
        max_mask = abs(votes - max_vote) < vote_tolerance;
        max_count = sum(max_mask);

        POS_consensus(obj_idx) = 0;
        BND_consensus(obj_idx) = 0;
        NEG_consensus(obj_idx) = 0;

        if max_count > 1
            BND_consensus(obj_idx) = 1;
        else
            [~,decision] = max(votes);

            if decision == 1
                POS_consensus(obj_idx) = 1;
            elseif decision == 2
                BND_consensus(obj_idx) = 1;
            else
                NEG_consensus(obj_idx) = 1;
            end
        end

        resolved(obj_idx) = true;
    end

    POS_T = POS_consensus;
    BND_T = BND_consensus;
    NEG_T = NEG_consensus;
end

%% 11. 安全兜底处理

if any(~resolved)

    unresolved_indices = find(~resolved);

    POS_T(unresolved_indices) = 0;
    BND_T(unresolved_indices) = 1;
    NEG_T(unresolved_indices) = 0;

    POS_consensus = POS_T;
    BND_consensus = BND_T;
    NEG_consensus = NEG_T;

    resolved(unresolved_indices) = true;
end

%% 12. 最终共识结果 Eq.(37)

POS_T = POS_consensus;
BND_T = BND_consensus;
NEG_T = NEG_consensus;

%% 13. 根据最终共识结果重新计算专家支持度

Su_final = calc_support(POS,BND,NEG,POS_T,BND_T,NEG_T,theta_expert,k_group);

%% 14. 最终综合得分 Eq.(40)

score_raw = nan(duixiang,1);

for i = 1:duixiang

    xbar_i = nan(1,k_group);

    for p = 1:k_group
        if ~isempty(S_ALL{p})
            xbar_i(p) = S_ALL{p}(i);
        end
    end

    if POS_T(i) == 1
        expert_mask = POS(:,i)' == 1;
    elseif BND_T(i) == 1
        expert_mask = BND(:,i)' == 1;
    elseif NEG_T(i) == 1
        expert_mask = NEG(:,i)' == 1;
    else
        expert_mask = false(1,k_group);
    end

    expert_mask = expert_mask & isfinite(xbar_i);

    if ~any(expert_mask)
        continue;
    end

    denominator = sum(Su_final(expert_mask));
    numerator = sum(Su_final(expert_mask) .* xbar_i(expert_mask));

    if denominator > eps
        score_raw(i) = numerator / denominator;
    else
        fallbackDen = sum(theta_expert(expert_mask));

        if fallbackDen > eps
            score_raw(i) = sum(theta_expert(expert_mask) .* xbar_i(expert_mask)) / fallbackDen;
        else
            score_raw(i) = mean(xbar_i(expert_mask));
        end
    end
end

%% 15. 最终排序

pos_idx = find(POS_T == 1);
bnd_idx = find(BND_T == 1);
neg_idx = find(NEG_T == 1);

pos_sorted = sort_indices_by_score(pos_idx,score_raw);
bnd_sorted = sort_indices_by_score(bnd_idx,score_raw);
neg_sorted = sort_indices_by_score(neg_idx,score_raw);

final_sorted_sequence = [pos_sorted(:);bnd_sorted(:);neg_sorted(:)];

%% 16. 计算ER、IER、PR、RR、F1和OD

C_true = find(POS0 == 1);
notC_true = find(NEG0 == 1);

POS_idx = find(POS_T == 1);
BND_idx = find(BND_T == 1);
NEG_idx = find(NEG_T == 1);

C_to_POS = numel(intersect(C_true,POS_idx));
C_to_NEG = numel(intersect(C_true,NEG_idx));
notC_to_POS = numel(intersect(notC_true,POS_idx));

ER = (C_to_NEG + notC_to_POS) / duixiang;

sigma = 0.8;

den_IER = numel(POS_idx) + numel(BND_idx);

if den_IER > 0
    IER = sigma * ER + (1-sigma) * numel(BND_idx) / den_IER;
else
    IER = sigma * ER;
end

den_PR = C_to_POS + notC_to_POS;

if den_PR > 0
    PR = C_to_POS / den_PR;
else
    PR = 0;
end

den_RR = C_to_NEG + C_to_POS;

if den_RR > 0
    RR = C_to_POS / den_RR;
else
    RR = 0;
end

if PR + RR > 0
    F1 = 2 * PR * RR / (PR + RR);
else
    F1 = 0;
end

rank_pos = zeros(1,duixiang);

for r = 1:numel(final_sorted_sequence)
    rank_pos(final_sorted_sequence(r)) = r;
end

J1 = intersect(C_true,NEG_idx);
J2 = intersect(notC_true,POS_idx);

if isempty(J1)
    avg_J1 = 0;
else
    avg_J1 = mean(rank_pos(J1));
end

if isempty(J2)
    avg_J2 = 0;
else
    avg_J2 = mean(rank_pos(J2));
end

if isempty(notC_true)
    o_JN = 0;
else
    o_JN = min(rank_pos(notC_true));
end

if isempty(C_true)
    o_JP = 0;
else
    o_JP = max(rank_pos(C_true));
end

OD = (numel(J1) + numel(J2)) / duixiang * (max(avg_J1,o_JN) - min(avg_J2,o_JP));

%% 17. 输出结果

POS_T
BND_T
NEG_T

Re = [sum(POS_T),sum(BND_T),sum(NEG_T)]

so = final_sorted_sequence'

W_num

fprintf('ER = %.6f\n',ER);
fprintf('IER = %.6f\n',IER);
fprintf('PR = %.6f\n',PR);
fprintf('RR = %.6f\n',RR);
fprintf('F1 = %.6f\n',F1);
fprintf('OD = %.6f\n',OD);

%% 18. 局部函数

function sortedIdx = sort_indices_by_score(indices,score_raw)

    if isempty(indices)
        sortedIdx = indices;
        return;
    end

    values = score_raw(indices);
    values(~isfinite(values)) = -Inf;

    [~,ord] = sort(values,'descend');
    sortedIdx = indices(ord);
end

function Su = calc_support(POS,BND,NEG,POS_T,BND_T,NEG_T,theta_expert,k_group)

    Su = zeros(1,k_group);

    for k = 1:k_group

        POS_k = POS(k,:) == 1;
        num_POS = sum(POS_k);

        if num_POS > 0
            sP = sum(POS_k & (POS_T == 1)) / num_POS;
        else
            sP = 0;
        end

        BND_k = BND(k,:) == 1;
        num_BND = sum(BND_k);

        if num_BND > 0
            sB = sum(BND_k & (BND_T == 1)) / num_BND;
        else
            sB = 0;
        end

        NEG_k = NEG(k,:) == 1;
        num_NEG = sum(NEG_k);

        if num_NEG > 0
            sN = sum(NEG_k & (NEG_T == 1)) / num_NEG;
        else
            sN = 0;
        end

        Su(k) = theta_expert(k) * (sP + sB + sN) / 3;
    end
end

function Influence = calc_influence_unresolved(POS,BND,NEG,theta_expert,resolved,valid_mask,k_group)

    Influence = zeros(1,k_group);
    tie_tol = 1e-12;

    for k = 1:k_group

        valid_k = valid_mask{k}(:);
        candidate_indices = find((~resolved) & valid_k);

        if isempty(candidate_indices)
            Influence(k) = 0;
            continue;
        end

        psi_sum = 0;

        for q = 1:numel(candidate_indices)

            i = candidate_indices(q);

            G_P = sum(theta_expert .* POS(:,i)');
            G_B = sum(theta_expert .* BND(:,i)');
            G_N = sum(theta_expert .* NEG(:,i)');

            G_vec = [G_P,G_B,G_N];
            G_i = max(G_vec);
            Gamma_i = abs(G_vec - G_i) < tie_tol;

            theta_minus = theta_expert;
            theta_minus(k) = 0;

            G_P_minus = sum(theta_minus .* POS(:,i)');
            G_B_minus = sum(theta_minus .* BND(:,i)');
            G_N_minus = sum(theta_minus .* NEG(:,i)');

            G_minus_vec = [G_P_minus,G_B_minus,G_N_minus];
            G_i_minus = max(G_minus_vec);
            Gamma_i_minus = abs(G_minus_vec - G_i_minus) < tie_tol;

            if any(Gamma_i_minus ~= Gamma_i)
                psi = 1;
            elseif G_i > eps
                psi = (G_i - G_i_minus) / G_i;
            else
                psi = 0;
            end

            psi_sum = psi_sum + psi;
        end

        Influence(k) = psi_sum / numel(candidate_indices);
    end
end

function [idx,L] = kmeans_mahal_global(X,k_class,opts)

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'useBiasCov')
        opts.useBiasCov = false;
    end

    if ~isfield(opts,'MaxIter')
        opts.MaxIter = 200;
    end

    if ~isfield(opts,'Replicates')
        opts.Replicates = 100;
    end

    if ~isfield(opts,'Start')
        opts.Start = 'plus';
    end

    if ~isfield(opts,'Seed')
        opts.Seed = 1;
    end

    if ~isfield(opts,'reg')
        opts.reg = 1e-6;
    end

    X = double(X);

    if any(~isfinite(X(:)))
        error('kmeans_mahal_global输入X中包含NaN或Inf。');
    end

    [n,p] = size(X);

    if n < k_class
        error('样本数量%d小于聚类数量%d。',n,k_class);
    end

    nUnique = size(unique(X,'rows'),1);

    if nUnique < k_class
        error('当前只有%d个不同数据点，但聚类数量为%d。',nUnique,k_class);
    end

    if opts.useBiasCov
        S = cov(X,1);
    else
        S = cov(X,0);
    end

    S = (S + S') / 2;

    [V,D] = eig(S);
    V = real(V);
    lambda = real(diag(D));
    lambda(~isfinite(lambda)) = 0;

    positiveLambda = lambda(lambda > 0);

    if isempty(positiveLambda)
        scaleCov = 1;
    else
        scaleCov = max(positiveLambda);
    end

    eigFloor = max(opts.reg * scaleCov,1e-10);

    lambdaReg = lambda;
    lambdaReg(lambdaReg < eigFloor) = eigFloor;

    Wwhite = V * diag(1 ./ sqrt(lambdaReg)) * V';

    mu = mean(X,1);
    Z = (X - mu) * Wwhite;

    rng(opts.Seed,'twister');

    idx = kmeans(Z,k_class,'Distance','sqeuclidean','Start',opts.Start,'Replicates',opts.Replicates,'MaxIter',opts.MaxIter,'EmptyAction','singleton','Display','off');

    L = nan(k_class,p);

    for c = 1:k_class

        mask = idx == c;

        if any(mask)
            L(c,:) = mean(X(mask,:),1);
        end
    end

    if any(~isfinite(L(:)))
        error('马氏距离K-means得到空聚类或无效聚类中心。');
    end
end
