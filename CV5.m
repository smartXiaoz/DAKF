close all;clear all; clc
N = 500;     %total number of time steps
n = 3;
m = 1;
T = 1;
F = [1 T T^2/2;
    0 1 T
    0 0 1]; % 模型
H = [0 1 0];
epsilon = 1e-6;
del_T=0.1;


xe_ini =ones(n,1)*5; %估计初始值
xe_MEE = xe_ini;
F_factor = 0.1;
Ke_est_MEE =0;Ke_est_MCC =0;Ke_est_MSE =0;Ke_est_MCC1 =0;
KeR_MEE = 0; KeR_MCC = 0; KeR_MSE = 0;KeR_MCC1 = 0;


Iter_num = 0;
sigma_MCC = 2; 
sigma_MEE= 1.5; 

sigma_MCC1 = 10; 
sigma_MEE1= 6; 
for mm = 1 : Iter_num
    Ke_est_MEE =0;Ke_est_MCC =0;Ke_est_MSE =0;Ke_est_MCC1 =0;
KeR_MEE = 0; KeR_MCC = 0; KeR_MSE = 0;KeR_MCC1 = 0;
    q = randn(n,N) * 0.1; % 过程噪声
    %% 观测噪声r
%     r1 = randn(m,N)*0.1; r2 = randn(m,N)*10;
%     r = randn(m,N);
%     for ii = 1:m
%         vp = rand(1,N);
%         r(ii,:) = (vp<=0.8).* r1(ii,:) +  (vp>0.8).*r2(ii,:);
%     end
    
    
    NN = 1000000;
%     NN = N;
    r = randn(m,NN) * 0.1;
%     mu1 = -5;
%     mu2 = 5;
    
    v1=randn(m,NN)*0.1 ; v2=randn(m,NN)*5 ;v3 = rand(m,NN)*2-1;v4 =  sin(randn(m,NN));
    rp=rand(1,NN);  
    v4=randn(m,NN)*40;
    %case1
%     r = v1;
    % case 2
%     r = (rp<=0.999).*v1 + (rp>0.999).*v4;
    % case 3
%     r = (rp<=0.8).*v1 + (rp>0.8).*v2;
%     v1 = randn(m,NN)*0.1+ 5;
%     v2 = randn(m,NN)*0.1 -5;
%      r = (rp<=0.5).*v1 + (rp>0.5).*v2;
%     % case 4
%     r = (rp<=0.8).*v1 + (rp>0.8).*v3;
    % case 5
 
    v1 = randn(m,NN)*0.1+ 0.5;
    v2 = randn(m,NN)*5 -0.5;
    v4 = gamrnd(4,5,m, NN);
    r = (rp<=0.8).*v1 + (rp>0.8 & rp<0.95).*v2 + (rp>0.95).*v4;
    %case5
%     v1 = randn(m,NN)*0.1+ 0.5;
%     v2 = randn(m,NN)*10 + 10;
%     v4 = gamrnd(0.3,5,m, NN)*0.1;
%     r = (rp<=0.8).*v1 + (rp>0.8 & rp<0.9).*v2 + (rp>0.9).*v4;

%      r = (rp<=0.5).*v1 + (rp>0.5).*v2;

%     r = (rp<=0.8).*v1 + (rp>0.8).*v4;
%     r = (rp<=0.99).*v1 + (rp>0.99).*v4;

%     r = raylrnd(1,m, NN);
%     r = gamrnd(0.3,5,m, NN)*0.1;
%     r = sin(randn(m,NN));
%     r = v1;
%     mu1 = 0;
%     mu2 = 0;
%     r = stblrnd(1.5,0,1.1,0,m,NN)*0.1;%alpha 分布噪声
%     r = stblrnd(1.2,0,3,0,m,NN);%alpha 分布噪声
    %% ***
    gmm=fitgmdist(r', 2,'RegularizationValue', 1e-8);
    
    R_pai = gmm.ComponentProportion;
    RR1 = gmm.Sigma(:,:,1);
    RR2 = gmm.Sigma(:,:,2);
    mu1 = gmm.mu(1);
    mu2 = gmm.mu(2);
   

    
    Q = q * q'/N;
    R  = r * r'/NN;
%     R  = diag(diag(r * r'/N));
    xx = zeros(n, N); %状态值
    yy = zeros(m, N); %观测值

    xx(:, 1) = [0 0 1];
    for ii = 2 : N
        xx(:, ii) = F * xx( :, ii-1 ) + q( :, ii-1 ); %true state update process
        yy(:, ii) = H * xx( :, ii ) + r( :, ii );%true measurement update process
    end
        
     %% MEE_KF循环计算    
    Pk_MEE = eye(n) * 1;
    xe_MEE = xe_ini;
    xx1 = xe_ini;
    Err_MEE_KF(mm,1) = ( xe_ini - xx( :, 1 ) )'*( xe_ini - xx( :, 1 ) );
    b_MEE(1) = 0;
%     tic
%     for uuu = 1:100
    for ii = 2:N      
     %% MEE定点迭代计算 xe
        yy1 = yy( :, ii );
        xx1 = xe_MEE( :, ii-1 );       
        [xeii_t,Pk_MEE,e,Phi_k,Psi_k,lambda,Ke] = MEE_KF_F(F,xx1,Pk_MEE,H,yy1,Q,R,sigma_MEE);       
        xe_MEE( :, ii ) = xeii_t;
        %% 计算误差值
        Err_MEE_KF(mm,ii) = (xe_MEE( :, ii ) - xx( :, ii ))'*(xe_MEE( :, ii ) - xx( :, ii ));       
    end
    
    for ii = 2:N      
     %% MEE定点迭代计算 xe
        yy1 = yy( :, ii );
        xx1 = xe_MEE( :, ii-1 );       
        [xeii_t,Pk_MEE,e,Phi_k,Psi_k,lambda,Ke] = MEE_KF_F(F,xx1,Pk_MEE,H,yy1,Q,R,sigma_MEE1);       
        xe_MEE( :, ii ) = xeii_t;
        %% 计算误差值
        Err_MEE_KF_1(mm,ii) = (xe_MEE( :, ii ) - xx( :, ii ))'*(xe_MEE( :, ii ) - xx( :, ii ));       
    end
%     end
%     toc
    %% MCC_KF循环计算
    xe_MCC = xe_ini;
    xx3 = xe_ini;
    Pk_MCC = eye(n) * 1;
%     Pk_MCC = diag([900 900 4 4]) ;
    Err_MCC_KF(mm,1) = ( xe_MCC( :,1 ) - xx( :, 1 ) )'*( xe_MCC( :,1 ) - xx( :, 1 ) );
    for ii = 2:N
        %% MCC定点迭代计算 xe3
        yy3 = yy( :, ii);      
        [xx3,Pk_MCC,b,C,Ke_MCC,M_MCC] = MCC_KF_WJX(F,xx3,Pk_MCC,H,yy3,Q,R,sigma_MCC);        
        xe_MCC( :, ii ) = xx3;
        %% 计算误差值
        Err_MCC_KF(mm,ii) = ( xe_MCC( :, ii ) - xx( :, ii ) )'*( xe_MCC( :, ii ) - xx( :, ii ) );        
    end
    
    %% MCC_KF循环计算
    xe_MCC = xe_ini;
    xx3 = xe_ini;
    Pk_MCC = eye(n) * 1;
%     Pk_MCC = diag([900 900 4 4]) ;
    Err_MCC_KF(mm,1) = ( xe_MCC( :,1 ) - xx( :, 1 ) )'*( xe_MCC( :,1 ) - xx( :, 1 ) );
    for ii = 2:N
        %% MCC定点迭代计算 xe3
        yy3 = yy( :, ii);      
        [xx3,Pk_MCC,b,C,Ke_MCC,M_MCC] = MCC_KF_WJX(F,xx3,Pk_MCC,H,yy3,Q,R,sigma_MCC1);        
        xe_MCC( :, ii ) = xx3;
        %% 计算误差值
        Err_MCC_KF_1(mm,ii) = ( xe_MCC( :, ii ) - xx( :, ii ) )'*( xe_MCC( :, ii ) - xx( :, ii ) );        
    end
    
    %% MSE_KF
    xe_MSE = xe_ini;    
    Pk2 = eye(n) * 1;
    Err_MSE_KF(mm,1) = ( xe_MSE - xx( :, 1 ) )'*( xe_MSE - xx( :, 1 ) );
    for ii = 2:N            
        Pke2 = F * Pk2 * F' + Q;                           
        G_MSE = Pke2 * H' * inv(H*Pke2*H' + R);
        xe_MSE( :, ii ) = F * xe_MSE( :, ii-1 ) + G_MSE*( yy( :, ii )-H*F*xe_MSE( :, ii-1 ) );
        Pk2 = Pke2-G_MSE*H*Pke2;        
        Err_MSE_KF(mm,ii) = (xe_MSE( :, ii )-xx( :, ii ))'*(xe_MSE( :, ii )-xx( :, ii ));
    end
%     tic
    %% EM-KF
%     for uuu = 1 : 100
    epsilo = 0.00001;
    xe_EM = xe_ini;
    Pkk_EM = eye(n)*1;       
    Err_EM_KF(mm,1) = ( xe_EM - xx( :, 1 ) )'*( xe_EM - xx( :, 1 ) );
    RRR = RR1*R_pai(1) + RR2*R_pai(2);
    for ii = 2:N
        yy1 = yy(:,ii);
        xke = F * xe_EM(:,ii-1);
        Pke = F * Pkk_EM * F' + Q;
        Bp = chol(Pke)';
        Br = chol(RRR)';
        BB = zeros(n+m);
        BB(1:n,1:n) = Bp;
        BB(n+1:n+m,n+1:n+m) = Br;
        DD = inv(BB) * [xke;yy1];
        WW = inv(BB) * [eye(n);H];        
        past = xke;
        xkk = xke;
        eey = yy1 - H * xke;
            ee = DD - WW * xkk;                      
            Cx = eye(n);                    
            Cy = eye(m);
            ee_test(:,ii) = ee(n+1:m+n);
            PP1 = H*Pke*H' + RR1;
            PP2 = H*Pke*H' + RR2;

            t1(ii) = 0;
            t2(ii) = 0;
            if ii > 30 
                vik1 = R_pai(1) * 1 / (2 * pi)^(m/2) / (abs(det(PP1)))^0.5 * exp(-0.5 * (eey-mu1)'*inv(PP1)*(eey-mu1));
                vik2 = R_pai(2) * 1 / (2 * pi)^(m/2) / (abs(det(PP2)))^0.5 * exp(-0.5 * (eey-mu2)'*inv(PP2)*(eey-mu2));      
                
                sumaa = (vik1 + vik2);
                if sumaa == 0
                    t1(ii) = 1;
                    t2(ii) = 0;
                else
                    t1(ii) = vik1   / sumaa;
                    t2(ii) = vik2   / sumaa;
                end
%                 tmprr1=sqrt(RR1)./sqrt(R);
%                 tmprr2=sqrt(RR2)./sqrt(R);
                Cy = ((t1(ii) * sqrt(RR2) + t2(ii) * sqrt(RR1)) ./ sqrt(RRR)); 
%                 Cy = t1(ii) ./ RR1 + t2(ii) ./ RR2;
                Cy = diag(diag(Cy));
            end
            R_hat = Br * inv(Cy) * Br';
            yy1 = yy1 - (t1(ii)*mu1 + t2(ii)*mu2);
%             Pke_hat = Bp * inv(Cx) * Bp';
            Gk_EM = Pke * H' *inv(H*Pke*H' + R_hat);
            xkk = xke + Gk_EM * (yy1 - H *  xke);
            xe_EM(:,ii) = xkk;
        tmp1 = (eye(n) - Gk_EM*H);
        Pkk_EM = tmp1*Pke*tmp1' + Gk_EM*R*Gk_EM';
        Err_EM_KF(mm,ii) = (xkk - xx(:,ii))'*(xkk - xx(:,ii));       
%         if mm == Iter_num
%             KeR_MCC = (1-F_factor)*KeR_MCC + F_factor* Gk_EM * r( :, ii )*r( :, ii )' * Gk_EM'; 
%             Ke_est_MCC = (1-F_factor)*Ke_est_MCC + F_factor*Gk_EM;  
%         end
    end
    
%     end
%     toc

%     Err_MCC_TH(mm) = fun_steady_state_error_KF(H,F,Q,Ke_est_MCC,KeR_MCC);
   
    
    disp(mm)
end
% Err_MCC_TH = fun_steady_state_error_KF(H,F,Q,Ke_est_MCC,KeR_MCC);
% Err_MCC_TH1 = fun_steady_state_error_KF(H,F,Q,Ke_est_MCC1,KeR_MCC1);
% figure, hold on;
% wid = 2.5;
% plot(10*log10(mean(Err_MSE_KF)),'m','LineWidth',wid);
% plot(10*log10(mean(Err_MCC_KF)),'c','LineWidth',wid);
% plot(10*log10(mean(Err_MEE_KF)),'r','LineWidth',wid);
% plot(10*log10(mean(Err_EM_KF)),'b','LineWidth',wid);
% 
% % plot(10*log10(mean(Err_MCC_TH)*ones(1,N)),'b--','LineWidth',wid);
% 
% h=legend('KF', 'MCKF', 'MEEKF', 'EMKF1','TH-EMKF1');
% set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
% set(gca,'fontsize',24);
% xlabel('Iterations','FontSize',24);
% ylabel('MSD(dB)','FontSize',24);

% 
Err_MSE_KF1 = sqrt(sum(Err_MSE_KF)/Iter_num);
Err_MCC_KF1 = sqrt(sum(Err_MCC_KF)/Iter_num);
Err_MCC_KF2 = sqrt(sum(Err_MCC_KF_1)/Iter_num);
Err_MEE_KF1 = sqrt(sum(Err_MEE_KF)/Iter_num);
Err_MEE_KF2 = sqrt(sum(Err_MEE_KF_1)/Iter_num);
Err_EM_KF1 = sqrt(sum(Err_EM_KF)/Iter_num);
% Err_MCC_TH1 = sqrt(sum(Err_MCC_TH)/Iter_num);

mark_idx = 1 : 100 : N;
figure, hold on;
box on
wid = 1;
plot(20*log10((Err_MSE_KF1)),'m','LineWidth',wid);
plot(20*log10((Err_MCC_KF1)),'c','LineWidth',wid);
plot(20*log10((Err_MEE_KF1)),'r','LineWidth',wid);
plot(20*log10((Err_EM_KF1)),'b-','LineWidth',wid);

% plot(20*log10(Err_MCC_TH1*ones(1,N)),'g.-','LineWidth',1,'MarkerIndices',mark_idx,'markersize',10);

h=legend('KF', 'MCKF \sigma=2', 'MEEKF \sigma=1.5', 'Proposed DAKF');
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal','orientation','horizontal');
set(gca,'fontsize',24);
xlabel('Iterations (case5)','FontSize',24);
ylabel('RMSE(dB)','FontSize',24);

A_MSE= 20*log10(mean(Err_MSE_KF1(:,200:end)))
A_MCC_2= 20*log10(mean(Err_MCC_KF1(:,200:end)))
A_MCC_10= 20*log10(mean(Err_MCC_KF2(:,200:end)))
A_MEE_15= 20*log10(mean(Err_MEE_KF1(:,200:end)))
A_MEE_6= 20*log10(mean(Err_MEE_KF2(:,200:end)))
A_EM = 20*log10(mean(Err_EM_KF1(:,200:end)))
% A_EMTH=20*log10(Err_MCC_TH1)
% A_MSE= (mean(Err_MSE_KF1(:,200:end)))
% A_MCC= (mean(Err_MCC_KF1(:,200:end)))
% A_MEE= (mean(Err_MEE_KF1(:,200:end)))
% A_EM = (mean(Err_EM_KF1(:,200:end)))
% A_EMTH=(Err_MCC_TH1)

% Err_MSE_KF2 = mean(sqrt(min(Err_MSE_KF(:,200:end))));
% Err_MCC_KF2 = min(mean(sqrt(Err_MCC_KF(:,200:end)),2));
% Err_MEE_KF2 = min(mean(sqrt(Err_MEE_KF(:,200:end)),2));
% Err_EM_KF2 = min(mean(sqrt(Err_EM_KF(:,200:end)),2));
% Err_MCC_TH2 = min(mean(sqrt(Err_MCC_TH),2));
% 
% 
% A1 = (mean(Err_MSE_KF1)) - Err_MSE_KF2
% A2 = (mean(Err_MCC_KF1)) - Err_MCC_KF2
% A3 = (mean(Err_MEE_KF1)) - Err_MEE_KF2
% A4 = (mean(Err_EM_KF1)) - Err_EM_KF2
% A5 = (mean(Err_MCC_TH1)) - Err_MCC_TH2