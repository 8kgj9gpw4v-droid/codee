clc;
clear;
close all;

%% 1. 基本参数

k_group  = 5;      % 专家组 / 数据组数量
duixiang = 20;     % 区域 / 备选对象数量
shuxing  = 5;      % 属性数量

rho = 0.1;

% 论文案例结果输出文件
%
% 仅输出论文案例正文 / Table V 实际展示的结果：
%   1) 专家覆盖系数 tau_k
%   2) 各专家属性权重
%   3) 各专家 POS/BND/NEG 数量
%   4) 最终 POS/BND/NEG 区域集合
%   5) 最终排序

outputResultFile = fullfile(pwd, 'CR_TWGDM_PaperCaseResults.xlsx');

%% 2. 输入文件

fileList = {
    'D:\Desktop\AllRegionAnalysisOutput\ValidRate_Output\D1_ValidRate_AllRegions.xlsx'
    'D:\Desktop\AllRegionAnalysisOutput\ValidRate_Output\D2_ValidRate_AllRegions.xlsx'
    'D:\Desktop\AllRegionAnalysisOutput\ValidRate_Output\D3_ValidRate_AllRegions.xlsx'
    'D:\Desktop\AllRegionAnalysisOutput\ValidRate_Output\D4_ValidRate_AllRegions.xlsx'
    'D:\Desktop\AllRegionAnalysisOutput\ValidRate_Output\D5_ValidRate_AllRegions.xlsx'
};

% 只有以下 5 个属性列参与计算。
% 第 6 列“decision”会被自动忽略。

colNames = {
    'range0_有效率%'
    'range1_有效率%'
    'range2_有效率%'
    '最大间距率%'
    '三圆交叉残差均值_m'
};

nFiles = numel(fileList);

if nFiles ~= k_group
    error('输入文件数量 (%d) 必须等于 k_group (%d)。', nFiles, k_group);
end

%% 3. 读取并预处理数据

input_data = cell(1, nFiles);
valid_mask = cell(1, nFiles);

for f = 1:nFiles

    if ~isfile(fileList{f})
        error('输入文件不存在：\n%s', fileList{f});
    end

    T = readtable(fileList{f}, 'VariableNamingRule', 'preserve');

    if height(T) ~= duixiang
        error('D%d 中共有 %d 行数据，但 duixiang = %d。\n程序默认 Excel 第 1 行对应区域 1，第 20 行对应区域 20。', f, height(T), duixiang);
    end

    missingCols = setdiff(colNames, T.Properties.VariableNames);

    if ~isempty(missingCols)
        error('D%d 缺少以下必要属性列：\n%s', f, strjoin(missingCols, ', '));
    end

    % 只读取参与计算的 5 个属性
    data = zeros(duixiang, shuxing);

    for j = 1:shuxing

        col = T.(colNames{j});

        % 将文本 / 字符串 / 单元格等数据转换为数值型

        col = column_to_double(col);

        % 属性 1-3：
        % 有效率统一转换到 [0,1]

        if j <= 3

            % 有效率缺失表示不存在有效样本，因此记为 0
            col(~isfinite(col)) = 0;

            % 同时兼容以下两种数据格式：
            %
            % 80, 75, 60
            %
            % 或：
            %
            % 0.80, 0.75, 0.60

            finiteVals = abs(col(isfinite(col)));

            if ~isempty(finiteVals) && max(finiteVals) > 1.5
                col = col / 100;
            end

            % 数值保护：将结果限制在 [0,1]
            col(col < 0) = 0;
            col(col > 1) = 1;

            data(:,j) = col;

            continue;
        end

        % 属性 4-5：
        % 成本型属性转换为效益型标准化值
        %
        % 只有当属性 1-3 全部为 0 时，
        % 才认为该区域为无效区域。
        %
        % 如果某区域本身有效，但成本型属性缺失，
        % 则其标准化后的效益值记为 0。

        first3ZeroMask = all(data(:,1:3) == 0, 2);

        finiteMask = ~first3ZeroMask & isfinite(col);

        finiteVals = col(finiteMask);

        % 默认处理：
        %
        % 有效区域 + 成本属性缺失 -> 0
        % 无效区域 -> 0

        newCol = zeros(size(col));

        if ~isempty(finiteVals)

            vmax = max(finiteVals);
            vmin = min(finiteVals);

            tolerance = eps(max(1, max(abs(finiteVals))));

            if abs(vmax - vmin) <= tolerance

                % 所有可用属性值完全相同时，
                % 将这些有效值统一标准化为 1

                newCol(finiteMask) = 1;

            else

                % 成本型属性 -> 效益型属性标准化：
                %
                % 最小成本 -> 1
                % 最大成本 -> 0

                newCol(finiteMask) = (vmax - col(finiteMask)) ./ (vmax - vmin);

            end
        end

        data(:,j) = newCol;
    end

    % 有效区域判定：
    % range0、range1、range2 中至少有一个不为 0

    valid_mask{f} = ~all(data(:,1:3) == 0, 2);

    % 最终数据安全检查

    Xcheck = data(valid_mask{f}, :);

    if any(~isfinite(Xcheck(:)))

        [rr, cc] = find(~isfinite(Xcheck), 1, 'first');

        validIdxTmp = find(valid_mask{f});

        error('D%d 经过数据预处理后仍然包含 NaN/Inf。\n区域 = %d，属性 = %s', f, validIdxTmp(rr), colNames{cc});
    end

    input_data{f} = data;
end

%% 4. 初始化变量

POS = zeros(k_group, duixiang);
BND = zeros(k_group, duixiang);
NEG = zeros(k_group, duixiang);

W_num     = nan(k_group, shuxing);
PP        = cell(1, k_group);
IDX_ALL   = cell(1, k_group);
S_ALL     = cell(1, k_group);
D_STATE   = false(k_group, duixiang);
C_DOM_ALL = cell(1, k_group);
ALPHA_ALL = cell(1, k_group);
BETA_ALL  = cell(1, k_group);

% 熵权法中间量

SUMX_ALL      = cell(1, k_group);
P_ENTROPY_ALL = cell(1, k_group);
H_ENTROPY_ALL = cell(1, k_group);
D_ENTROPY_ALL = cell(1, k_group);

% K-means / 状态集关键中间量

CENTER_ALL       = cell(1, k_group);
CENTER_SCORE_ALL = cell(1, k_group);
STATE_CLASS_ALL  = nan(1, k_group);
KMEANS_LOSS_ALL  = nan(1, k_group);
KMEANS_INFO_ALL  = cell(1, k_group);

%% 5. 对每一个专家 / 数据集分别进行三支决策分类

for p = 1:k_group

    Xfull = input_data{p};

    validRowMask = valid_mask{p};

    validIdx = find(validRowMask);

    nValid = numel(validIdx);

    if nValid < 2
        error('D%d 仅包含 %d 个有效区域，至少需要 2 个有效区域。', p, nValid);
    end

    Pro_full   = nan(duixiang,1);
    idx_full   = nan(duixiang,1);
    xbar_full  = nan(duixiang,1);
    alpha_full = nan(duixiang,1);
    beta_full  = nan(duixiang,1);

    X = Xfull(validRowMask, 1:shuxing);

    if any(~isfinite(X(:)))
        error('D%d：进行马氏距离 K-means 聚类前，X 中包含 NaN 或 Inf。', p);
    end

    %% 5.1 熵权法计算属性权重 Eq.(8)-Eq.(10)

    sumX = sum(X, 1);

    % Eq.(8) 的分母必须为正。
    % 论文没有给出“整列和为 0”时的替代定义，因此这里不再静默修改分母。

    if any(sumX <= 0)
        badCols = find(sumX <= 0);
        error('D%d：Eq.(8) 无法计算，属性列 %s 的列和 <= 0。', p, mat2str(badCols));
    end

    % Eq.(8)

    P = X ./ sumX;

    % Eq.(9)：严格实现 p*ln(p)，并按论文规定 0*ln(0)=0。
    % 不使用 log(P+eps)，也不再对 H 做 [0,1] 截断。

    plnp = zeros(size(P));
    positiveMaskP = (P > 0);
    plnp(positiveMaskP) = P(positiveMaskP) .* log(P(positiveMaskP));

    entropy_k = 1 / log(nValid);

    H = -entropy_k * sum(plnp, 1);

    % Eq.(10)：严格使用 d_j = 1-H_j。
    % 不再使用 max(1-H,0) 截断。

    d = 1 - H;

    weightDen = sum(d);

    % 论文 Eq.(10) 在分母为 0 时没有给出退化规则，
    % 因此这里直接报错，而不是自动退化为等权重。

    if abs(weightDen) <= eps
        error('D%d：Eq.(10) 权重分母 sum(1-H_j) 为 0，论文未定义该退化情形。', p);
    end

    W = d / weightDen;

    W_num(p,:) = W;

    % 保存熵权法中间量，方便与论文案例逐项核对

    P_full_entropy = nan(duixiang, shuxing);
    P_full_entropy(validRowMask,:) = P;

    SUMX_ALL{p}      = sumX;
    P_ENTROPY_ALL{p} = P_full_entropy;
    H_ENTROPY_ALL{p} = H;
    D_ENTROPY_ALL{p} = d;

    %% 5.2 加权相对优势关系 Eq.(11)-Eq.(15)

    C_dom = zeros(nValid, nValid);

    for i = 1:nValid

        for r = 1:nValid

            diff_ir = X(i,:) - X(r,:);

            % Eq.(11)：对象 i 相对于对象 r 的加权优势

            A_ir = sum(W .* max(diff_ir, 0));

            % Eq.(12)：对象 i 相对于对象 r 的加权劣势

            L_ir = sum(W .* max(-diff_ir, 0));

            H_ir = A_ir + L_ir;

            % Eq.(13)：计算相对优势度

            if H_ir > eps

                C_dom(i,r) = A_ir / H_ir;

            else

                % 两个对象完全一致时，
                % 相对优势度定义为 0.5

                C_dom(i,r) = 0.5;

            end
        end
    end

    C_DOM_ALL{p} = C_dom;

    %% 5.3 稳健全局马氏距离 K-means 聚类 Eq.(16)

    k_class = 2;

    options = struct();

    options.useBiasCov = false;
    options.MaxIter    = 200;
    options.Replicates = 1000;
    options.Start      = 'plus';
    options.Seed       = 1;
    options.reg        = 1e-6;

    [idx, L, cluster_loss, cluster_info] = kmeans_mahal_global(X, k_class, options);

    %% 5.4 状态集合 D^k Eq.(17)

    % 根据聚类中心的加权综合得分，
    % 找出表现更优的聚类作为状态集合 D^k

    center_score = L * W';

    [~, state_class] = max(center_score);

    D_valid = (idx == state_class);

    D_STATE(p, validIdx(D_valid)) = true;

    CENTER_ALL{p}       = L;
    CENTER_SCORE_ALL{p} = center_score;
    STATE_CLASS_ALL(p)  = state_class;
    KMEANS_LOSS_ALL(p)  = cluster_loss;
    KMEANS_INFO_ALL{p}  = cluster_info;

    %% 5.5 条件概率 Eq.(18)

    Pro_valid = zeros(nValid,1);

    for i = 1:nValid

        % 找出相对于对象 i 具有优势关系的对象

        dominance_class = find(C_dom(:,i) >= 0.5);

        if isempty(dominance_class)

            Pro_valid(i) = 0;

        else

            % 计算优势类中属于状态集合 D^k 的对象比例

            Pro_valid(i) = sum(D_valid(dominance_class)) / numel(dominance_class);

        end
    end

    % 将有效区域的结果填回原始 20 个区域的位置

    Pro_full(validRowMask) = Pro_valid;
    idx_full(validRowMask) = idx;

    PP{p}      = Pro_full';
    IDX_ALL{p} = idx_full;

    %% 5.6 相对损失函数及决策阈值 Eq.(19)

    % 每个区域的属性加权综合值

    xbar = X * W';

    xbar_max = max(xbar);
    xbar_min = min(xbar);

    alpha = zeros(nValid,1);
    beta  = zeros(nValid,1);

    for i = 1:nValid

        % 与最佳对象之间的距离

        A = xbar_max - xbar(i);

        % 与最差对象之间的距离

        B = xbar(i) - xbar_min;

        den_alpha = (1-rho)*A + rho*B;

        den_beta = rho*A + (1-rho)*B;

        if den_alpha > eps

            alpha(i) = ((1-rho)*A) / den_alpha;

        else

            alpha(i) = 0;

        end

        if den_beta > eps

            beta(i) = (rho*A) / den_beta;

        else

            beta(i) = 0;

        end
    end

    xbar_full(validRowMask) = xbar;

    alpha_full(validRowMask) = alpha;
    beta_full(validRowMask)  = beta;

    S_ALL{p}     = xbar_full;
    ALPHA_ALL{p} = alpha_full;
    BETA_ALL{p}  = beta_full;

    %% 5.7 三支决策 Eq.(20)
    %
    % 条件概率 >= alpha       -> POS（正域）
    % beta < 条件概率 < alpha -> BND（边界域）
    % 条件概率 <= beta        -> NEG（负域）

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

%% 6. 专家信任系数 Eq.(22)
%
% 这里使用每个专家有效区域数占全部区域数的比例
% 作为专家信任系数

theta_expert = zeros(1, k_group);

for p = 1:k_group

    theta_expert(p) = sum(valid_mask{p}) / duixiang;

end

%% 7. 初始共识集合 Eq.(23)

% 根据专家信任系数对各专家的 POS/BND/NEG 判断进行加权

weighted_POS0 = theta_expert * POS;

weighted_BND0 = theta_expert * BND;

weighted_NEG0 = theta_expert * NEG;

consensus_tol = 1e-12;

POS_T = zeros(1, duixiang);
BND_T = zeros(1, duixiang);
NEG_T = zeros(1, duixiang);

%% POS 初始共识

max_POS0 = max(weighted_POS0);

if max_POS0 > consensus_tol

    % 将 POS 加权支持度达到最大值的对象
    % 作为初始 POS 共识对象

    POS_T(abs(weighted_POS0 - max_POS0) < consensus_tol) = 1;

end

%% BND 初始共识

max_BND0 = max(weighted_BND0);

if max_BND0 > consensus_tol

    % 将 BND 加权支持度达到最大值的对象
    % 作为初始 BND 共识对象

    BND_T(abs(weighted_BND0 - max_BND0) < consensus_tol) = 1;

end

%% NEG 初始共识

max_NEG0 = max(weighted_NEG0);

if max_NEG0 > consensus_tol

    % 将 NEG 加权支持度达到最大值的对象
    % 作为初始 NEG 共识对象

    NEG_T(abs(weighted_NEG0 - max_NEG0) < consensus_tol) = 1;

end

POS_consensus = POS_T;
BND_consensus = BND_T;
NEG_consensus = NEG_T;

% 保存 Eq.(23) 得到的初始一致性集合，供论文案例核对

POS_INITIAL = POS_T;
BND_INITIAL = BND_T;
NEG_INITIAL = NEG_T;

%% 8. 迭代冲突消解 Eq.(25)-Eq.(36)

Su        = zeros(1, k_group);
Influence = zeros(1, k_group);
Do        = zeros(1, k_group);

Co = nan(duixiang, 1);

% 保存每轮 Eq.(25)-Eq.(36) 的关键中间量

ITER_LOG = struct('iteration', {}, 'unresolved', {}, 'Su', {}, 'Influence', {}, 'Do', {}, 'WP', {}, 'WB', {}, 'WN', {}, 'PhiP', {}, 'PhiB', {}, 'PhiN', {}, 'Conflict', {}, 'target_objs', {}, 'target_decision', {}, 'POS_consensus', {}, 'BND_consensus', {}, 'NEG_consensus', {});

iteration = 0;

% resolved(i)=true 表示第 i 个对象已经完成冲突消解

resolved = false(duixiang, 1);

% 理论上每次至少解决一个对象，因此最多进行 duixiang+1 次

max_iterations = duixiang + 1;

while any(~resolved) && iteration < max_iterations

    iteration = iteration + 1;

    unresolved_indices = find(~resolved);

    if isempty(unresolved_indices)
        break;
    end

    % Eq.(25)-Eq.(26)：计算专家支持度

    Su = calc_support(POS, BND, NEG, POS_T, BND_T, NEG_T, theta_expert, k_group);

    % Eq.(27)-Eq.(30)：计算专家对未解决对象的影响度

    Influence = calc_influence_unresolved(POS, BND, NEG, theta_expert, resolved, valid_mask, k_group);

    % Eq.(31)：计算专家重要度
    %
    % 专家重要度 = 专家支持度 × 专家影响度

    Do = Su .* Influence;

    % Eq.(32)-Eq.(34)：计算每个未解决对象的冲突度

    Co = nan(duixiang,1);

    WP_iter   = nan(duixiang,1);
    WB_iter   = nan(duixiang,1);
    WN_iter   = nan(duixiang,1);
    PhiP_iter = nan(duixiang,1);
    PhiB_iter = nan(duixiang,1);
    PhiN_iter = nan(duixiang,1);

    for i = unresolved_indices'

        % POS 的加权投票值

        weighted_P = sum(Do .* POS(:,i)');

        % BND 的加权投票值

        weighted_B = sum(Do .* BND(:,i)');

        % NEG 的加权投票值

        weighted_N = sum(Do .* NEG(:,i)');

        W_region = [weighted_P, weighted_B, weighted_N];

        total_weight = sum(W_region);

        if total_weight <= eps

            % 如果所有有效投票权重都为 0，
            % 则认为 POS/BND/NEG 三种决策完全不确定，
            % 即概率均为 1/3

            phi = [1/3, 1/3, 1/3];

        else

            % 归一化得到三种决策的比例

            phi = W_region / total_weight;

        end

        WP_iter(i)   = weighted_P;
        WB_iter(i)   = weighted_B;
        WN_iter(i)   = weighted_N;
        PhiP_iter(i) = phi(1);
        PhiB_iter(i) = phi(2);
        PhiN_iter(i) = phi(3);

        phi_nonzero = phi(phi > 0);

        % 使用归一化信息熵表示冲突程度
        % Co 越小，说明决策冲突越小

        Co(i) = -sum(phi_nonzero .* log(phi_nonzero)) / log(3);
    end

    % Eq.(35)：选择当前冲突度最小的对象

    Co_unresolved = Co(unresolved_indices);

    if any(~isfinite(Co_unresolved))
        error('第 %d 次迭代过程中出现非有限冲突度。', iteration);
    end

    min_co_val = min(Co_unresolved);

    co_tolerance = 1e-10;

    % 如果多个对象的最小冲突度相同，
    % 则这些对象在本轮中同时处理

    target_objs = unresolved_indices(abs(Co_unresolved - min_co_val) < co_tolerance);

    if isempty(target_objs)
        error('冲突消解过程中未能选择待处理对象。');
    end

    % Eq.(36)：对选中的对象重新进行三支决策分类

    target_decision = strings(numel(target_objs),1);

    for t = 1:numel(target_objs)

        obj_idx = target_objs(t);

        % 基于专家重要度计算 POS/BND/NEG 三类加权投票

        vote_POS = sum(Do .* POS(:,obj_idx)');

        vote_BND = sum(Do .* BND(:,obj_idx)');

        vote_NEG = sum(Do .* NEG(:,obj_idx)');

        votes = [vote_POS, vote_BND, vote_NEG];

        max_vote = max(votes);

        vote_tolerance = 1e-10;

        max_mask = abs(votes - max_vote) < vote_tolerance;

        max_count = sum(max_mask);

        % 先清除当前对象原有的共识分类

        POS_consensus(obj_idx) = 0;
        BND_consensus(obj_idx) = 0;
        NEG_consensus(obj_idx) = 0;

        % 如果最大投票结果不唯一，即出现并列，
        % 则将该对象划分到 BND（边界域）

        if max_count > 1

            BND_consensus(obj_idx) = 1;
            target_decision(t) = "BND";

        else

            [~, decision] = max(votes);

            if decision == 1

                POS_consensus(obj_idx) = 1;
                target_decision(t) = "POS";

            elseif decision == 2

                BND_consensus(obj_idx) = 1;
                target_decision(t) = "BND";

            else

                NEG_consensus(obj_idx) = 1;
                target_decision(t) = "NEG";

            end
        end

        % 标记该对象已经完成冲突消解

        resolved(obj_idx) = true;
    end

    % 保存本轮所有关键中间量

    ITER_LOG(iteration).iteration       = iteration;
    ITER_LOG(iteration).unresolved      = unresolved_indices(:);
    ITER_LOG(iteration).Su              = Su;
    ITER_LOG(iteration).Influence       = Influence;
    ITER_LOG(iteration).Do              = Do;
    ITER_LOG(iteration).WP              = WP_iter;
    ITER_LOG(iteration).WB              = WB_iter;
    ITER_LOG(iteration).WN              = WN_iter;
    ITER_LOG(iteration).PhiP            = PhiP_iter;
    ITER_LOG(iteration).PhiB            = PhiB_iter;
    ITER_LOG(iteration).PhiN            = PhiN_iter;
    ITER_LOG(iteration).Conflict        = Co;
    ITER_LOG(iteration).target_objs     = target_objs(:);
    ITER_LOG(iteration).target_decision = target_decision;
    ITER_LOG(iteration).POS_consensus   = POS_consensus;
    ITER_LOG(iteration).BND_consensus   = BND_consensus;
    ITER_LOG(iteration).NEG_consensus   = NEG_consensus;

    % 更新当前共识结果，
    % 用于下一轮计算专家支持度

    POS_T = POS_consensus;
    BND_T = BND_consensus;
    NEG_T = NEG_consensus;
end

%% 9. 安全兜底处理
%
% 如果达到最大迭代次数后仍存在未解决对象，
% 则统一将其划分到 BND（边界域），
% 表示其决策仍然具有较强不确定性。

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

%% 10. 最终共识结果 Eq.(37)

POS_T = POS_consensus;
BND_T = BND_consensus;
NEG_T = NEG_consensus;

%% 11. 根据最终共识结果重新计算专家支持度

Su_final = calc_support(POS, BND, NEG, POS_T, BND_T, NEG_T, theta_expert, k_group);

%% 12. 最终综合得分 Eq.(40)

score_raw = nan(duixiang,1);

for i = 1:duixiang

    xbar_i = nan(1, k_group);

    % 读取每个专家对区域 i 的属性综合得分

    for p = 1:k_group

        if ~isempty(S_ALL{p})

            xbar_i(p) = S_ALL{p}(i);

        end
    end

    % 选择与最终三支决策结果一致的专家
    %
    % 如果最终结果为 POS：
    % 只使用把该区域判定为 POS 的专家。
    %
    % 如果最终结果为 BND：
    % 只使用把该区域判定为 BND 的专家。
    %
    % 如果最终结果为 NEG：
    % 只使用把该区域判定为 NEG 的专家。

    if POS_T(i) == 1

        expert_mask = (POS(:,i)' == 1);

    elseif BND_T(i) == 1

        expert_mask = (BND(:,i)' == 1);

    elseif NEG_T(i) == 1

        expert_mask = (NEG(:,i)' == 1);

    else

        expert_mask = false(1, k_group);

    end

    % 同时要求该专家对当前区域具有有效综合得分

    expert_mask = expert_mask & isfinite(xbar_i);

    if ~any(expert_mask)
        continue;
    end

    % 优先使用最终专家支持度作为加权系数

    denominator = sum(Su_final(expert_mask));

    numerator = sum(Su_final(expert_mask) .* xbar_i(expert_mask));

    if denominator > eps

        score_raw(i) = numerator / denominator;

    else

        % 备用方案 1：
        % 如果专家支持度之和为 0，
        % 则改用专家信任系数 theta_expert 进行加权

        fallbackDen = sum(theta_expert(expert_mask));

        if fallbackDen > eps

            score_raw(i) = sum(theta_expert(expert_mask) .* xbar_i(expert_mask)) / fallbackDen;

        else

            % 备用方案 2：
            % 如果专家信任系数之和也为 0，
            % 则直接计算这些专家综合得分的算术平均值

            score_raw(i) = mean(xbar_i(expert_mask));

        end
    end
end

%% 13. 最终排序
%
% 第一层排序优先级：
%
% POS > BND > NEG
%
% 即：
% 正域对象排在最前，
% 边界域对象排在中间，
% 负域对象排在最后。
%
% 在相同三支决策类别内部：
%
% 按 FinalScore / score_raw 从大到小排序。

pos_idx = find(POS_T == 1);

bnd_idx = find(BND_T == 1);

neg_idx = find(NEG_T == 1);

% POS 内部按综合得分降序排列

pos_sorted = sort_indices_by_score(pos_idx, score_raw);

% BND 内部按综合得分降序排列

bnd_sorted = sort_indices_by_score(bnd_idx, score_raw);

% NEG 内部按综合得分降序排列

neg_sorted = sort_indices_by_score(neg_idx, score_raw);

% 最终顺序：
% POS -> BND -> NEG

final_sorted_sequence = [
    pos_sorted(:)
    bnd_sorted(:)
    neg_sorted(:)
];

%% 14. 论文案例展示量输出
%
% 只输出论文案例正文与 Table V 实际展示的量：
%   1) Eq.(22) 专家覆盖系数 tau_k
%   2) Eq.(10) 各专家属性权重
%   3) 各专家 POS/BND/NEG 数量
%   4) Eq.(37) 最终 POS/BND/NEG 区域集合
%   5) 最终排序

%% 14.1 各专家三支决策数量

expert_partition_counts = [sum(POS,2), sum(BND,2), sum(NEG,2)];

%% 14.2 最终三支决策区域集合

final_POS_regions = find(POS_T == 1);

final_BND_regions = find(BND_T == 1);

final_NEG_regions = find(NEG_T == 1);

%% 14.3 命令窗口输出

fprintf('\n============================================================\n');
fprintf('CR-TWGDM 论文案例结果\n');
fprintf('============================================================\n');

fprintf('\n1. 专家覆盖系数 tau_k：\n');
disp(theta_expert);

fprintf('2. 各专家属性权重 W：\n');
disp(W_num);

fprintf('3. 各专家三支决策数量 [POS  BND  NEG]：\n');
disp(expert_partition_counts);

fprintf('4. 最终三支决策区域：\n');

fprintf('POS = {');
fprintf('u%d ', final_POS_regions);
fprintf('}\n');

fprintf('BND = {');
fprintf('u%d ', final_BND_regions);
fprintf('}\n');

fprintf('NEG = {');
fprintf('u%d ', final_NEG_regions);
fprintf('}\n');

fprintf('5. 最终排序：\n');

for r = 1:numel(final_sorted_sequence)

    if r < numel(final_sorted_sequence)
        fprintf('u%d > ', final_sorted_sequence(r));
    else
        fprintf('u%d\n', final_sorted_sequence(r));
    end
end

fprintf('============================================================\n');

%% 14.4 保存为论文成稿用 Excel

if isfile(outputResultFile)
    delete(outputResultFile);
end

% Sheet 1: 专家覆盖系数

T_tau = table((1:k_group)', theta_expert', 'VariableNames', {'Expert', 'Tau'});

writetable(T_tau, outputResultFile, 'Sheet', 'Coverage');

% Sheet 2: 各专家属性权重

T_weight = table((1:k_group)', W_num(:,1), W_num(:,2), W_num(:,3), W_num(:,4), W_num(:,5), 'VariableNames', {'Expert', 'a1', 'a2', 'a3', 'a4', 'a5'});

writetable(T_weight, outputResultFile, 'Sheet', 'Weights');

% Sheet 3: 各专家 POS/BND/NEG 数量

T_partition = table((1:k_group)', expert_partition_counts(:,1), expert_partition_counts(:,2), expert_partition_counts(:,3), 'VariableNames', {'Expert', 'POS', 'BND', 'NEG'});

writetable(T_partition, outputResultFile, 'Sheet', 'ExpertPartitions');

% Sheet 4: 最终分类

maxClassLen = max([numel(final_POS_regions), numel(final_BND_regions), numel(final_NEG_regions)]);

POS_col = strings(maxClassLen,1);
BND_col = strings(maxClassLen,1);
NEG_col = strings(maxClassLen,1);

for q = 1:numel(final_POS_regions)
    POS_col(q) = "u" + string(final_POS_regions(q));
end

for q = 1:numel(final_BND_regions)
    BND_col(q) = "u" + string(final_BND_regions(q));
end

for q = 1:numel(final_NEG_regions)
    NEG_col(q) = "u" + string(final_NEG_regions(q));
end

T_final_class = table(POS_col, BND_col, NEG_col, 'VariableNames', {'POS', 'BND', 'NEG'});

writetable(T_final_class, outputResultFile, 'Sheet', 'FinalClassification');

% Sheet 5: 最终排序

Rank = (1:numel(final_sorted_sequence))';
Region = "u" + string(final_sorted_sequence(:));

T_ranking = table(Rank, Region, 'VariableNames', {'Rank', 'Region'});

%% 局部函数：
% 将表格中的某一列统一转换为 double 类型

function x = column_to_double(col)

    if isnumeric(col) || islogical(col)

        % 数值型或逻辑型直接转换

        x = double(col);

    elseif iscell(col) || isstring(col) || ischar(col) || iscategorical(col)

        % 其他文本类数据统一转换为 string

        s = string(col);

        % 删除首尾空格

        s = strtrim(s);

        % 删除百分号

        s = erase(s, '%');

        % 判断常见的空值表示方式

        emptyMask = (s == "") | ismissing(s) | strcmpi(s, "nan") | strcmpi(s, "na") | strcmpi(s, "n/a") | strcmpi(s, "null") | strcmpi(s, "none");

        % 统一设置为缺失值

        s(emptyMask) = missing;

        % 转换为 double

        x = str2double(s);

    else

        error('不支持的列数据类型：%s', class(col));

    end

    % 强制转换成列向量

    x = x(:);
end

%% 局部函数：
% 根据综合得分从高到低对区域编号进行排序
% NaN / Inf 等非有限值统一排到最后

function sortedIdx = sort_indices_by_score(indices, score_raw)

    if isempty(indices)

        sortedIdx = indices;
        return;

    end

    values = score_raw(indices);

    % 非有限得分视为 -Inf，
    % 从而保证这些区域排在最后

    values(~isfinite(values)) = -Inf;

    [~, ord] = sort(values, 'descend');

    sortedIdx = indices(ord);
end

%% Eq.(25)-Eq.(26)：
% 专家支持度计算

function Su = calc_support(POS, BND, NEG, POS_T, BND_T, NEG_T, theta_expert, k_group)

    Su = zeros(1, k_group);

    for k = 1:k_group

        %% POS 支持度

        % 专家 k 判断为 POS 的对象

        POS_k = (POS(k,:) == 1);

        num_POS = sum(POS_k);

        if num_POS > 0

            % 专家 k 的 POS 判断中，
            % 与当前共识 POS 一致的比例

            sP = sum(POS_k & (POS_T == 1)) / num_POS;

        else

            sP = 0;

        end

        %% BND 支持度

        % 专家 k 判断为 BND 的对象

        BND_k = (BND(k,:) == 1);

        num_BND = sum(BND_k);

        if num_BND > 0

            % 专家 k 的 BND 判断中，
            % 与当前共识 BND 一致的比例

            sB = sum(BND_k & (BND_T == 1)) / num_BND;

        else

            sB = 0;

        end

        %% NEG 支持度

        % 专家 k 判断为 NEG 的对象

        NEG_k = (NEG(k,:) == 1);

        num_NEG = sum(NEG_k);

        if num_NEG > 0

            % 专家 k 的 NEG 判断中，
            % 与当前共识 NEG 一致的比例

            sN = sum(NEG_k & (NEG_T == 1)) / num_NEG;

        else

            sN = 0;

        end

        %% 专家总体支持度

        % 将 POS、BND、NEG 三类支持度取平均，
        % 再乘以专家自身信任系数

        Su(k) = theta_expert(k) * (sP + sB + sN) / 3;
    end
end

%% Eq.(27)-Eq.(30)：
% 专家对当前未解决对象的影响度

function Influence = calc_influence_unresolved(POS, BND, NEG, theta_expert, resolved, valid_mask, k_group)

    Influence = zeros(1, k_group);

    tie_tol = 1e-12;

    for k = 1:k_group

        % 专家 k 对哪些区域具有有效数据

        valid_k = valid_mask{k}(:);

        % 只考虑：
        % 1. 当前尚未解决；
        % 2. 对专家 k 而言有效
        % 的区域

        candidate_indices = find((~resolved) & valid_k);

        if isempty(candidate_indices)

            Influence(k) = 0;
            continue;

        end

        psi_sum = 0;

        for i = candidate_indices'

            % Eq.(27)：
            % 保留专家 k 时的群体加权决策结果

            G_P = sum(theta_expert .* POS(:,i)');

            G_B = sum(theta_expert .* BND(:,i)');

            G_N = sum(theta_expert .* NEG(:,i)');

            G_vec = [G_P, G_B, G_N];

            % 当前最大的群体加权支持值

            G_i = max(G_vec);

            % 判断哪些类别与最大支持值并列

            Gamma_i = abs(G_vec - G_i) < tie_tol;

            % Eq.(28)：
            % 去掉专家 k 后重新计算群体加权决策结果

            theta_minus = theta_expert;

            % 将专家 k 的权重置为 0，
            % 相当于从群体决策中移除该专家

            theta_minus(k) = 0;

            G_P_minus = sum(theta_minus .* POS(:,i)');

            G_B_minus = sum(theta_minus .* BND(:,i)');

            G_N_minus = sum(theta_minus .* NEG(:,i)');

            G_minus_vec = [
                G_P_minus
                G_B_minus
                G_N_minus
            ]';

            G_i_minus = max(G_minus_vec);

            Gamma_i_minus = abs(G_minus_vec - G_i_minus) < tie_tol;

            % Eq.(29)：
            % 计算专家 k 对对象 i 的局部影响程度

            % 如果删除专家 k 后，
            % 最大决策类别集合发生变化，
            % 说明该专家对最终决策具有关键影响，
            % 因此 psi = 1

            if any(Gamma_i_minus ~= Gamma_i)

                psi = 1;

            elseif G_i > eps

                % 如果最大类别没有变化，
                % 则使用最大支持值下降比例衡量影响程度

                psi = (G_i - G_i_minus) / G_i;

            else

                psi = 0;

            end

            psi_sum = psi_sum + psi;
        end

        % Eq.(30)：
        % 对所有候选对象的局部影响程度取平均，
        % 得到专家 k 的总体影响度

        Influence(k) = psi_sum / numel(candidate_indices);
    end
end

%% 稳健全局马氏距离 K-means 聚类函数
%
% 基本思想：
%
% 1. 根据全部样本计算全局协方差矩阵；
% 2. 对协方差矩阵进行特征值正则化；
% 3. 使用协方差矩阵对白化空间进行变换；
% 4. 在白化后的空间中执行普通欧氏距离 K-means；
% 5. 该过程等价于在原始空间中使用全局马氏距离聚类。

function [idx, L, loss, info] = kmeans_mahal_global(X, k_class, opts)

    % 默认参数设置

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'useBiasCov')
        opts.useBiasCov = false;
    end

    if ~isfield(opts, 'MaxIter')
        opts.MaxIter = 200;
    end

    if ~isfield(opts, 'Replicates')
        opts.Replicates = 100;
    end

    if ~isfield(opts, 'Start')
        opts.Start = 'plus';
    end

    if ~isfield(opts, 'Seed')
        opts.Seed = 1;
    end

    if ~isfield(opts, 'reg')
        opts.reg = 1e-6;
    end

    % 基本数据检查

    X = double(X);

    if any(~isfinite(X(:)))

        error('kmeans_mahal_global：X 中包含 NaN 或 Inf。');

    end

    [n, p] = size(X);

    if n < k_class

        error('样本数量 (%d) 小于聚类数量 k_class (%d)。', n, k_class);

    end

    % 计算不同样本点的数量

    nUnique = size(unique(X, 'rows'), 1);

    if nUnique < k_class

        error('当前只有 %d 个不同的数据点，但 k_class = %d。', nUnique, k_class);

    end

    % 计算全局协方差矩阵

    if opts.useBiasCov

        % 使用有偏协方差估计

        S = cov(X, 1);

    else

        % 使用无偏协方差估计

        S = cov(X, 0);

    end

    % 强制保证矩阵对称

    S = (S + S') / 2;

    % 特征值正则化
    %
    % 防止协方差矩阵奇异或接近奇异，
    % 从而提高马氏距离计算的数值稳定性

    [V, D] = eig(S);

    V = real(V);

    lambda = real(diag(D));

    % 非有限特征值统一设为 0

    lambda(~isfinite(lambda)) = 0;

    positiveLambda = lambda(lambda > 0);

    if isempty(positiveLambda)

        scaleCov = 1;

    else

        scaleCov = max(positiveLambda);

    end

    % 设置最小允许特征值

    eigFloor = max(opts.reg * scaleCov, 1e-10);

    lambdaReg = lambda;

    % 小于阈值的特征值统一提升到 eigFloor

    lambdaReg(lambdaReg < eigFloor) = eigFloor;

    % 重构正则化后的协方差矩阵

    Sreg = V * diag(lambdaReg) * V';

    Sreg = (Sreg + Sreg') / 2;

    % 白化变换
    %
    % 将原始数据映射到一个新的空间，
    % 使得该空间中的欧氏距离对应原始空间中的
    % 全局马氏距离。

    Wwhite = V * diag(1 ./ sqrt(lambdaReg)) * V';

    % 原始数据均值

    mu = mean(X,1);

    % 对数据中心化后执行白化变换

    Z = (X - mu) * Wwhite;

    % 在白化空间中执行 K-means

    % 固定随机种子，
    % 保证程序每次运行时具有可重复性

    rng(opts.Seed, 'twister');

    [idx, ~, sumd] = kmeans(Z, k_class, 'Distance', 'sqeuclidean', 'Start', opts.Start, 'Replicates', opts.Replicates, 'MaxIter', opts.MaxIter, 'EmptyAction', 'singleton', 'Display', 'off');

    % 在原始属性空间中重新计算各聚类中心

    L = nan(k_class, p);

    for c = 1:k_class

        mask = (idx == c);

        if any(mask)

            % 聚类中心为该类所有原始样本的属性均值

            L(c,:) = mean(X(mask,:),1);

        end
    end

    if any(~isfinite(L(:)))

        error('马氏距离 K-means 得到了空聚类或无效的聚类中心。');

    end

    % 总聚类损失

    loss = sum(sumd);

    % 保存诊断信息
    %
    % 这些信息不参与后续主要决策，
    % 主要用于检查聚类和协方差正则化过程。

    info = struct();

    % 原始数据均值

    info.mean = mu;

    % 原始协方差矩阵

    info.originalCovariance = S;

    % 正则化后的协方差矩阵

    info.regularizedCovariance = Sreg;

    % 原始特征值

    info.originalEigenvalues = lambda;

    % 正则化后的特征值

    info.regularizedEigenvalues = lambdaReg;

    % 特征值最小阈值

    info.eigenvalueFloor = eigFloor;

    % 白化矩阵

    info.whiteningMatrix = Wwhite;

    % 聚类损失

    info.loss = loss;

    % 样本数量

    info.numberOfSamples = n;

    % 属性数量

    info.numberOfAttributes = p;

    % 不重复样本数量

    info.numberOfUniqueSamples = nUnique;
end
