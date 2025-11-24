function CCo = Kai2Opt3D1(FcS, phaseShift, G)
% varargin: FcS1aT ... FcS12aT

% === 准备工作 ===
% FcS = cat(3, varargin{:}); % [H, W, 12]
numPhase = numel(phaseShift) + 1; % +1 for phaseShift1=0
nHarm = 6; % 1~6 阶谐波

% === 构造相位矩阵 ===
phaseAll = [0; phaseShift(:)]; % 13×1
% P = exp(-1i * phaseAll); % 基波
% P(1) = 1; % phaseShift1 = 0 对应 exp(0)=1

% === 构造 M 矩阵 ===
M = zeros(numPhase - 1, 2 * nHarm); % 12×12
for k = 1:nHarm
    base = exp(-1i * k * phaseAll(2:end)) - exp(-1i * phaseAll(1));
    M(:, 2*k-1) = base;
    M(:, 2*k)   = conj(base);
end
Minv = inv(M); % 12×12

% === 合并 FcS 信号 ===
H = size(FcS, 1);
W = size(FcS, 2);
FcS = reshape(FcS, [], 12); % (H*W)×12

% 计算所有 12 通道调制结果（FiSMap1~6, FiSMam1~6）
Fi = FcS * Minv.'; % (H*W)×12
Fi = reshape(Fi, H, W, 12); % [H, W, 12]

% === 计算 CCo ===


% 奇偶索引：奇数=FiSMapN，偶数=FiSMamN
CCo = 0;
for k = 1:nHarm
    map = Fi(:, :, 2*k-1);
    mam = Fi(:, :, 2*k);
    CCo = CCo + abs(sum(sum(map .* conj(mam) .* G)));
end
end
