close all;clear all; clc
N = 5000;     %total number of time steps
n = 2;
m = 1;
F = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)]; % 模型
H = ones(m,n);
epsilon = 1e-6;
del_T=0.1;


xe_ini =ones(n,1)*4; %估计初始值
xe_MEE = xe_ini;
F_factor = 0.15;
Ke_est_MEE =0;Ke_est_MCC =0;Ke_est_MSE =0;Ke_est_MCC1 =0;
KeR_MEE = 0; KeR_MCC = 0; KeR_MSE = 0;KeR_MCC1 = 0;


Iter_num = 10;
sigma_MCC = 0.5; 
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
    
    
    NN = 10000;
    r = randn(m,NN) * 0.1;
%     mu1 = -5;
%     mu2 = 5;
    mu1 = 0;
    mu2 = 0;
    v1=randn(m,NN)*0.1 + mu2; v2=randn(m,NN)*10 + mu1;
    rp=rand(1,NN);  
    r = (rp<=0.8).*v1 + (rp>0.8).*v2;
    r = raylrnd(0.1,m, NN);
%     mu1 = 0;
%     mu2 = 0;
%     r = stblrnd(1.5,0,3,0,m,N)*0.5;%alpha 分布噪声
%     v = stblrnd(1.2,0,3,0,m,N);%alpha 分布噪声
    %% ***
    gmm=fitgmdist(r', 2);
    
    R_pai = gmm.ComponentProportion;
    RR1 = gmm.Sigma(:,:,1);
    RR2 = gmm.Sigma(:,:,2);
%     [R_pai,RR1,RR2] = EM_GMM2(r,mu1,mu2);
%     [R_pai,RR1,RR2,RR3] = EM_GMM3(r,mu1,mu2);
%     RR1 = diag(diag(RR1));
%     RR2 = diag(diag(RR2));
    
    Q = q * q'/N;
    R  = r * r'/NN;
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
        eey = yy1 - H * xke;
            ee = DD - WW * xkk;                      
            Cx = eye(n);                    
            Cy = eye(m);
            ee_test(:,ii) = ee(n+1:m+n);
            PP1 = H*Pke*H' + RR1;
            PP2 = H*Pke*H' + RR2;           

            if ii > 1 
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
                tmprr1=sqrt(RR1)./sqrt(R);
                tmprr2=sqrt(RR2)./sqrt(R);
                Cy = ((t1(ii) * sqrt(RR2) + t2(ii) * sqrt(RR1)) ./ sqrt(R)); 
%                 Cy = t1(ii) ./ RR1 + t2(ii) ./ RR2;
                Cy = diag(diag(Cy));
            end
            R_hat = Br * inv(Cy) * Br';
%             yy1 = yy1 - (t1(ii)*mu1 + t2(ii)*mu2);
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

    
    %% EM-KF
    epsilo = 0.00001;
    xe_EM = xe_ini;
    Pkk_EM = eye(n)*1;       
    Err_EM_KF1(mm,1) = ( xe_EM - xx( :, 1 ) )'*( xe_EM - xx( :, 1 ) );
    for ii = 2:N
        yy1 = yy(:,ii);
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
        eey = yy1 - H * xke;
            ee = DD - WW * xkk;                      
            Cx = eye(n);                    
            Cy = eye(m);
            ee_test(:,ii) = ee(n+1:m+n);
            PP1 = H*Pke*H' + RR1;
            PP2 = H*Pke*H' + RR2;           
            if ii > 1
               for jj = 1 : m
                   if abs(eey(jj)) > 10 * abs(PP1(jj,jj)).^0.5            
                        t2(ii) = 1;
                        t1(ii) = 0;
                    else
                         t2(ii) = 0;
                         t1(ii) = 1;
                    end
                Cy(jj,jj) = ((t1(ii) * sqrt(RR2(jj,jj)) + t2(ii) * sqrt(RR1(jj,jj))) ./ sqrt(R(jj,jj))); 
               end
            end
            R_hat = Br * inv(Cy) * Br';
%             yy1 = yy1 - (t1(ii)*mu1 + t2(ii)*mu2);
            Pke_hat = Bp * inv(Cx) * Bp';
            Gk_EM = Pke_hat * H' *inv(H*Pke_hat*H' + R_hat);
            xkk = xke + Gk_EM * (yy1 - H *  xke);
            xe_EM(:,ii) = xkk;
        tmp1 = (eye(n) - Gk_EM*H);
        Pkk_EM = tmp1*Pke*tmp1' + Gk_EM*R*Gk_EM';
        Err_EM_KF1(mm,ii) = (xkk - xx(:,ii))'*(xkk - xx(:,ii));       
        if mm == Iter_num
            F_factor = 0.12;
            KeR_MCC1 = (1-F_factor)*KeR_MCC1 + F_factor* Gk_EM * r( :, ii )*r( :, ii )' * Gk_EM'; 
            Ke_est_MCC1 = (1-F_factor)*Ke_est_MCC1 + F_factor*Gk_EM;  
        end
    end
    
    disp(mm)
end
Err_MCC_TH = fun_steady_state_error_KF(H,F,Q,Ke_est_MCC,KeR_MCC);
Err_MCC_TH1 = fun_steady_state_error_KF(H,F,Q,Ke_est_MCC1,KeR_MCC1);
figure, hold on;
wid = 2.5;
plot(10*log10(mean(Err_MSE_KF)),'m','LineWidth',wid);
plot(10*log10(mean(Err_MCC_KF)),'c','LineWidth',wid);
plot(10*log10(mean(Err_MEE_KF)),'r','LineWidth',wid);
plot(10*log10(mean(Err_EM_KF)),'b','LineWidth',wid);
plot(10*log10(mean(Err_EM_KF1)),'g','LineWidth',wid);
plot(10*log10(Err_MCC_TH*ones(1,N)),'b--','LineWidth',wid);
plot(10*log10(Err_MCC_TH1*ones(1,N)),'g--','LineWidth',wid);
h=legend('KF', 'MCKF', 'MEEKF', 'EMKF1','EMKF2','TH-EMKF1','TH-EMKF2');
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
set(gca,'fontsize',24);
xlabel('Iterations','FontSize',24);
ylabel('MSD(dB)','FontSize',24);
