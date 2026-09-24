function result = approximation(C, varargin)
% APPROXIMATION
% No normalization is performed.
%
% The first seven columns are used as conditional attributes by default.
% The last column, if present, is treated as the decision attribute and
% does not participate in approximation calculation.

    %% Read data

    if nargin == 0 || isempty(C)

        names = { ...
            'appendicitis处理.xls', ...
            'appendicitis处理.xlsx', ...
            'appendicitis.xls', ...
            'appendicitis.xlsx'};

        found = false;

        for k = 1:numel(names)
            if exist(names{k}, 'file') == 2
                C = readmatrix(names{k});
                found = true;
                break;
            end
        end

        if ~found
            error('未提供输入数据，且当前目录未找到 Appendicitis 数据文件。');
        end
    end

    validateattributes(C, ...
        {'numeric', 'logical'}, ...
        {'2d', 'nonempty'}, ...
        mfilename, 'C');

    C = double(C);
    raw = C;

    [m, totalColumns] = size(C);

    %% Parameters

    p = inputParser;

    addParameter(p, 'conditionColumns', [], ...
        @(x) isempty(x) || ...
        (isnumeric(x) && isvector(x) && all(x >= 1)));

    addParameter(p, 'alpha', 0.05, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

    addParameter(p, 'beta', 0.05, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

    addParameter(p, 'lambda', 0.50, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);

    addParameter(p, 'radius', 0.20, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    addParameter(p, 'eta', 0.20, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    addParameter(p, 'dengScale', [], ...
        @(x) isempty(x) || ...
        (isnumeric(x) && isscalar(x) && x > 0));

    parse(p, varargin{:});
    q = p.Results;

    %% Select conditional attributes

    if isempty(q.conditionColumns)

        if totalColumns == 8
            conditionColumns = 1:7;
        else
            conditionColumns = 1:totalColumns;
        end

    else
        conditionColumns = q.conditionColumns(:).';

        if any(conditionColumns > totalColumns)
            error('conditionColumns 超出输入矩阵的列数。');
        end
    end

    Xraw = C(:, conditionColumns);

    [m, pAttr] = size(Xraw);

    %% No normalization

    % 直接使用原始条件属性值
    X = Xraw;

    % Replace NaN and Inf values by the median of the corresponding column
    for a = 1:pAttr

        z = X(:, a);
        validMask = isfinite(z);

        if ~any(validMask)
            X(:, a) = 0;
        else
            z(~validMask) = median(z(validMask));
            X(:, a) = z;
        end
    end

    %% Notice

    if any(X(:) < 0 | X(:) > 1)
        warning(['当前未进行归一化，部分属性值不在 [0,1] 范围内。', ...
                 '若将原始属性值直接作为模糊隶属度，所得结果可能不再属于 [0,1]。']);
    end

    %% Initialize output matrices

    WenUpper = zeros(m, pAttr);
    WenLower = zeros(m, pAttr);

    YuanUpper = zeros(m, pAttr);
    YuanLower = zeros(m, pAttr);

    LDFRSUpper = zeros(m, pAttr);
    LDFRSLower = zeros(m, pAttr);

    DengUpper = zeros(m, pAttr);
    DengLower = zeros(m, pAttr);

    %% 1. Wen et al.: attribute-wise construction

    for a = 1:pAttr

        % Raw attribute values are used as membership degrees
        mu = X(:, a);

        D = abs(mu - mu.');

        scale = max(D(:));

        if scale <= eps
            scale = 1;
        end

        R = exp(-D / scale);

        % Upper approximation
        WenUpper(:, a) = ...
            min(1, max(min(R, mu.'), [], 2) + q.beta);

        % Lower approximation
        WenLower(:, a) = ...
            max(0, min(max(1 - R, mu.'), [], 2) - q.alpha);
    end

    %% 2. Yuan et al.: attribute-wise fuzzy neighborhood construction

    for a = 1:pAttr

        mu = X(:, a);

        D = abs(mu - mu.');

        N = max(0, 1 - D / q.radius);

        % Remove weak neighborhood relations
        N(N < q.lambda) = 0;

        % Upper approximation
        YuanUpper(:, a) = max(min(N, mu.'), [], 2);

        % Lower approximation
        YuanLower(:, a) = min(max(1 - N, mu.'), [], 2);
    end

    %% 3. Yu and Yao: attribute-wise logical-distance construction

    for a = 1:pAttr

        mu = X(:, a);

        D = abs(mu - mu.');

        logicalDistance = min(1, D / q.eta);

        % Upper approximation
        LDFRSUpper(:, a) = ...
            max(min(1, logicalDistance + mu.'), [], 2);

        % Lower approximation
        LDFRSLower(:, a) = ...
            min(max(1 - logicalDistance, mu.'), [], 2);
    end

    %% 4. Deng and Wang: joint construction based on all attributes

    % Calculate Euclidean distance using all conditional attributes
    squaredDistance = zeros(m, m);

    for a = 1:pAttr

        difference = X(:, a) - X(:, a).';

        squaredDistance = ...
            squaredDistance + difference .^ 2;
    end

    EuclideanDistance = sqrt(squaredDistance);

    % Scale parameter in R_H
    if isempty(q.dengScale)
        dengScale = pAttr;
    else
        dengScale = q.dengScale;
    end

    RH = 1 - EuclideanDistance / dengScale;

    % Keep the fuzzy relation within [0,1]
    RH = min(max(RH, 0), 1);

    % The same joint relation is used for every attribute
    for a = 1:pAttr

        mu = X(:, a);

        % Upper approximation
        DengUpper(:, a) = max(min(RH, mu.'), [], 2);

        % Lower approximation
        DengLower(:, a) = min(max(1 - RH, mu.'), [], 2);
    end

    %% Return results

    result = { ...
        raw, ...
        struct( ...
            'paper', 'Wen et al. (VWMFRS)', ...
            'upper', WenUpper, ...
            'lower', WenLower), ...
        struct( ...
            'paper', 'Yuan et al. (FNRS)', ...
            'upper', YuanUpper, ...
            'lower', YuanLower), ...
        struct( ...
            'paper', 'Yu and Yao (LDFRS)', ...
            'upper', LDFRSUpper, ...
            'lower', LDFRSLower), ...
        struct( ...
            'paper', 'Deng and Wang et al. (Euclidean-distance FRS)', ...
            'upper', DengUpper, ...
            'lower', DengLower) ...
        };
end