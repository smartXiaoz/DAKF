function [P_k,P_k_k_1,X_k]=kalman(F,T,H,Q,R,Z,X0,P0)
X_k_k_1 = F * X0;  %k-1时刻对k时刻x值的预测
P_k_k_1 = F*P0*F' + T*Q*T'; %k-1时刻对k时刻p值的预测
K_k = P_k_k_1 * H' * inv(H*P_k_k_1*H' + R);%k时刻kalman滤波增益
X_k = X_k_k_1+K_k*(Z - H*X_k_k_1);
P_k = P_k_k_1 - P_k_k_1* H' * inv(H*P_k_k_1*H' + R) * H * P_k_k_1;
end
