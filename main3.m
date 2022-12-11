close all;clear all; clc
N = 2000;     %total number of time steps
n = 2;
m = 1;
F = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)]; % 模型
H = ones(m,n);
epsilon = 1e-6;
del_T=0.1;


xe_ini =ones(n,1)*4; %估计初始值
xe_MEE = xe_ini;
F_factor = 0.15;
Ke_est_MEE =0;Ke_est_MCC =0;Ke_est_MSE =0;
KeR_MEE = 0; KeR_MCC = 0; KeR_MSE = 0;


Iter_num = 10;
sigma_MCC = 1; 
sigma_MEE=  0.5;

for mm = 1 : Iter_num
    q = randn(n,N) * 0.05; % 过程噪声
    %% 观测噪声r
%     r1 = randn(m,N)*0.1; r2 = randn(m,N)*10;
%     r = randn(m,N);
%     for ii = 1:m
%         vp = rand(1,N);
%         r(ii,:) = (vp<=0.8).* r1(ii,:) +  (vp>0.8).*r2(ii,:);
%     end
    
    
    
    r = randn(m,N) * 0.1;
%     mu1 = -5;
%     mu2 = 5;
    mu1 = 0;
    mu2 = 0;
    v1=randn(m,N)*0.1 + mu2; v2=randn(m,N)*20 + mu1;
    rp=rand(1,N);  
    r = (rp<=0.8).*v1 + (rp>0.8).*v2;
%     mu1 = 0;
%     mu2 = 0;
%     r = stblrnd(1.5,0,3,0,m,N)*0.5;%alpha 分布噪声
%     v = stblrnd(1.2,0,3,0,m,N);%alpha 分布噪声
    %% ***
    
        gmm=fitgmdist(r', 2);
    
    R_pai = gmm.ComponentProportion;
    RR1 = gmm.Sigma(:,:,1);
    RR2 = gmm.Sigma(:,:,2);
%     [R_pai,RR1,RR2,RR3] = EM_GMM3(r,mu1,mu2);
%     RR1 = diag(diag(RR1));
%     RR2 = diag(diag(RR2));
    
    Q = q * q'/N;
    R  = r * r'/N;
%     R  = diag(diag(r * r'/N));
    xx = zeros(n, N); %状态值
    yy = zeros(m, N); %观测值

    


    


    
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
    tic
    for ii = 2:N      
        %% MEE定点迭代计算 xe
        yy1 = yy( :, ii );
        xx1 = xe_MEE( :, ii-1 );       
        [xeii_t,Pk_MEE,e,Phi_k,Psi_k,lambda,Ke] = MEE_KF_F(F,xx1,Pk_MEE,H,yy1,Q,R,sigma_MEE);       
        xe_MEE( :, ii ) = xeii_t;
        %% 计算误差值
        Err_MEE_KF(mm,ii) = (xe_MEE( :, ii ) - xx( :, ii ))'*(xe_MEE( :, ii ) - xx( :, ii ));       
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

    %% EM-KF
    epsilo = 0.00001;
    xe_EM = xe_ini;
    Pkk_EM = eye(n)*1;  
       
    Err_EM_KF(mm,1) = ( xe_EM - xx( :, 1 ) )'*( xe_EM - xx( :, 1 ) );
    for ii = 2:N
        yy1 = yy(:,ii);
%         yy1 = yy1 - R_pai(1)*mu1 - R_pai(2)*mu2;
        xke = F * xe_EM(:,ii-1);
        Pke = F * Pkk_EM * F' + Q;
        Bp = chol(Pke)';
        Br = chol(R)';
        BB = zeros(n+m);
        BB(1:n,1:n) = Bp;
        BB(n+1:n+m,n+1:n+m) = Br;
        DD = inv(BB) * [xke;yy1];
        WW = inv(BB) * [eye(n);H];
        
        past = xke;
        xkk = xke;
%         eey = yy1 - H * xke;
%         while 1
            ee = DD - WW * xkk;                      
            Cx = eye(n);
            
            
            
            Cy = eye(m);
            ee_test(:,ii) = ee(n+1:m+n);
            eee = [xke;yy1]-[eye(n);H]*xkk;
            eex = ee(1:n);
            eey = eee(n+1:m+n);
%             eey = r(:,ii);
            eee_dd(:,ii) = r(:,ii) - eee(n+1:m+n);
            PP1 = H*Pke*H' + RR1;
            PP2 = H*Pke*H' + RR2;
            
            PPR = H*Pke*H' + R;
            PD1 = inv(Br)*RR1*inv(Br)';
            PD2 = inv(Br)*RR2*inv(Br)';
            if ii > 1
                for  i = 1 : m 
                    sigma = PP1(i,i);
                    Vik(i,1) = R_pai(1) * exp(-1 * (eey(i) - 0)^2/(2*sigma)) / sqrt(2*pi*sigma); 
                    sigma = PP2(i,i);    
                    Vik(i,2) = R_pai(2) * exp(-1 * (eey(i) - 0)^2/(2*sigma)) / sqrt(2*pi*sigma); 
                    sumaa = (Vik(i,1) + Vik(i,2));
                    
                    
                    if sumaa == 0
                        t1(i,ii) = 1;
                        t2(i,ii) = 0;
                    else
                        t1(i,ii) = Vik(i,1) / sumaa;
                        t2(i,ii) = Vik(i,2) / sumaa;
                    end
                    
                   
                    pp1 =  PP1(i,i) / eey(i)^2;
                    pp2 =  PP2(i,i) / eey(i)^2;
%                     t2(i,ii) = 1 - t1(i,ii);
%                     tmp1 = t1(i,ii) * (RR1(i,i));
%                     tmp2 = t2(i,ii) * (RR2(i,i));
%                     Cy(i,i) = (tmp1 + tmp2);
                    Cy(i,i) = ((t1(i,ii) * sqrt(RR2(i,i)) + t2(i,ii) * sqrt(RR1(i,i)))); 
                    Cy(i,i) = ((t1(i,ii) * sqrt(RR2(i,i)) + t2(i,ii) * sqrt(RR1(i,i)))/sqrt(R(i,i))); 
%                     Cy(i,i) = ((t1(i,ii) / sqrt(RR1(i,i)) + t2(i,ii) / sqrt(RR2(i,i)))*(R(i,i))); 
%                     Cy(i,i) = Vik(i,2) + Vik(i,1);
%                     Cy(i,i) = t1(i,ii) * sqrt(RR2(i,i)) + t2(i,ii) * sqrt(RR1(i,i));
%                     Cy(i,i) = ((t1(i,ii) * (RR2(i,i)) + t2(i,ii) * (RR1(i,i))))/R(i,i); 
%                     Cy(i,i) = (t1(i,ii) * sqrt(PP2(i,i)) + t2(i,ii) * sqrt(PP1(i,i)))/sqrt(R_pai(1)*PP1(i,i) + R_pai(2)*PP2(i,i));
                end  
            end

%             if ii > 1 
%                 vik1 = R_pai(1) * 1 / (2 * pi)^(m/2) / (abs(det(PP1)))^0.5 * exp(-0.5 * (eey-mu1)'*inv(PP1)*(eey-mu1));
%                 vik2 = R_pai(2) * 1 / (2 * pi)^(m/2) / (abs(det(PP2)))^0.5 * exp(-0.5 * (eey-mu2)'*inv(PP2)*(eey-mu2));
%                
%                 sumaa = (vik1 + vik2);
% %                 sumaa = (vik1*tmp1 + vik2*tmp2);
%                 if sumaa == 0
%                     t1(ii) = 1;
%                     t2(ii) = 0;
%                 else
%                     t1(ii) = vik1   / sumaa;
%                     t2(ii) = vik2   / sumaa;
% 
%                 end
%                 
% %                 if abs(sum(eey)) > 3 * sum(diag(PP2.^0.5))
% %                     t1(ii) = 1;
% %                     t2(ii) = 0;
% %                 else
% %                      t1(ii) = 0;
% %                      t2(ii) = 1;
% %                 end
%                 Cy = ((t1(ii) * sqrt(RR2) + t2(ii) * sqrt(RR1)) ./ sqrt(R)); 
% %                 PA1 = eye(m) * 5;
% %                 PA2 = eye(m) * 0.01;
% %                 Cy = ((t1(ii) * (PA2) + t2(ii) * (PA1))); 
% %                 Cy = ((t1(ii) * (PD2) + t2(ii) * (PD1)))./ (inv(Br)*R*inv(Br)');
% %                 Cy = (t1(ii) ./ (RR2) + t2(ii) ./ (RR2)) .* R;
%                 Cy = diag(diag(Cy));
%             end

            R_hat = Br * inv(Cy) * Br';


            yy1 = yy1 - (t1(ii)*mu1 + t2(ii)*mu2);
%             yy1 = yy1 - (R_pai(1)*mu1 + R_pai(2)*mu2);
            Pke_hat = Bp * inv(Cx) * Bp';
            Gk_EM = Pke_hat * H' *inv(H*Pke_hat*H' + R_hat);
            xkk = xke + Gk_EM * (yy1 - H *  xke);
            xe_EM(:,ii) = xkk;
        tmp1 = (eye(n) - Gk_EM*H);
        Pkk_EM = tmp1*Pke*tmp1' + Gk_EM*R*Gk_EM';
        Err_EM_KF(mm,ii) = (xkk - xx(:,ii))'*(xkk - xx(:,ii));
        
        if mm == Iter_num
            KeR_MCC = (1-F_factor)*KeR_MCC + F_factor* Gk_EM * r( :, ii )*r( :, ii )' * Gk_EM'; 
            Ke_est_MCC = (1-F_factor)*Ke_est_MCC + F_factor*Gk_EM;  
        end
    end

    disp(mm)
end
Err_MCC_TH = fun_steady_state_error_KF(H,F,Q,Ke_est_MCC,KeR_MCC);
figure, hold on;
plot(10*log10(mean(Err_MSE_KF)),'r');
plot(10*log10(mean(Err_MEE_KF)),'b');
plot(10*log10(mean(Err_EM_KF)),'k');
plot(10*log10(Err_MCC_TH*ones(1,N)),'g--');
xlabel('Iterations');ylabel('MSD(dB)');
