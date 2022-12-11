function [alphaK, miuK, sigmaK] = F_EM(noiseSeq, orders)

%% 输入 噪声序列noiseSeq， 混合阶数 orders
%% 输出 混合比例alphaK 相应地均值 miuK 方差 sigmaK 

miuK=[0 1];%一开始随机给的均值
sigmaK=[1 2];%一开始随机给的方差
alphaK=ones(1,orders)/orders;%一开始随机给的权重

LEN=length(noiseSeq);%信号长度

for mm = 1:10
    for ii = 1 : LEN
        vv = noiseSeq(ii);
        for kk = 1 : orders
            miu = miuK(kk);
            sigma = sigmaK(kk);
            Vik(ii, kk) = alphaK(kk) * exp(-1 * (vv - miu)^2/(2*sigma)) / sqrt(2*pi*sigma);            
        end
        Vik(ii, :) = Vik(ii, :) / (sum(Vik(ii, :)) + eps);        
    end
    NK = 0;
    for ii = 1 : LEN
        NK = NK + Vik(ii, :);
    end
    
    for kk = 1 : orders
        tmp1 = 0; 
        tmp2 = 0;
        for ii = 1 : LEN
            tmp1 = tmp1 + Vik(ii, kk) * noiseSeq(ii) / NK(kk);           
        end
        miuK(kk) = tmp1; % 更新第k个的均值        
        for ii = 1 : LEN
            tmp2 = tmp2 + Vik(ii, kk) * (noiseSeq(ii) - miuK(kk))^2 / NK(kk);           
        end
        sigmaK(kk) = tmp2;
        alphaK(kk) = NK(kk) / LEN;
    end
       
end