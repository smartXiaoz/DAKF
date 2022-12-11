close all;clear all; clc
N = 2000;     %total number of time steps
n = 2;
m = 2;
F = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)]; % 模型
H = ones(m,n);
epsilon = 1e-6;
del_T=0.1;


xe_ini =ones(n,1)*4; %估计初始值
xe_MEE = xe_ini;




Iter_num = 10;
sigma_MCC = 1; 
sigma_MEE=  0.2;

for mm = 1 : Iter_num
    q = randn(n,N) * 0.01; % 过程噪声
    %% 观测噪声r
%     r1 = randn(m,N)*0.1; r2 = randn(m,N)*10;
%     r = randn(m,N);
%     for ii = 1:m
%         vp = rand(1,N);
%         r(ii,:) = (vp<=0.8).* r1(ii,:) +  (vp>0.8).*r2(ii,:);
%     end
    
    
    
    r = randn(m,N) * 0.1;
    v1=randn(m,N)*0.1; v2=randn(m,N)*10;
    rp=rand(1,N);  
    r = (rp<=0.8).*v1 + (rp>0.8).*v2;
%     r = stblrnd(1.5,0,0.1,0,m,N) * 0.1;%alpha 分布噪声
    %% ***
    Q = q * q'/N;
    R  = r * r'/N;
    R  = diag(diag(r * r'/N));
    xx = zeros(n, N); %状态值
    yy = zeros(m, N); %观测值
    %% ***
    br = chol(R)';
    bbr = inv(br);
    
    for ir = 1 : N    
        rr(:,ir) = r(:,ir);
%         rr(:,ir) = bbr * r(:,ir);
    end
    for ir = 1:size(r,1)
        [alphaK(ir,:), miuK(ir,:), sigmaK(ir,:)] = F_EM(rr(ir,:),2);
    end   
    R_pai = mean(alphaK);    
    RR1 = diag(sigmaK(:,1));
    RR2 = (diag(sigmaK(:,2))); 

%     for ir = 1:size(r,1)
%         [alphaK(ir,:), miuK(ir,:), sigmaK(ir,:)] = F_EM(q(ir,:),2);
%     end   
%     Q_pai = mean(alphaK);    
%     QQ1 = diag(sigmaK(:,1));
%     QQ2 = diag(sigmaK(:,2)); 
    
    for ii = 2 : N
        xx(:, ii) = F * xx( :, ii-1 ) + q( :, ii-1 ); %true state update process
        yy(:, ii) = H * xx( :, ii ) + r( :, ii );%true measurement update process
    end
        tic
     %% MEE_KF循环计算    
    Pk_MEE = eye(n) * 1;
    xe_MEE = xe_ini;
    xx1 = xe_ini;
    Err_MEE_KF(mm,1) = ( xe_ini - xx( :, 1 ) )'*( xe_ini - xx( :, 1 ) );
    b_MEE(1) = 0;
    for ii = 2:N      
        %% MEE定点迭代计算 xe
        yy1 = yy( :, ii );
        xx1 = xe_MEE( :, ii-1 );       
        [xeii_t,Pk_MEE,e,Phi_k,Psi_k,lambda,Ke] = MEE_KF_F(F,xx1,Pk_MEE,H,yy1,Q,R,sigma_MEE);       
        xe_MEE( :, ii ) = xeii_t;
        %% 计算误差值
        Err_MEE_KF(mm,ii) = (xe_MEE( :, ii ) - xx( :, ii ))'*(xe_MEE( :, ii ) - xx( :, ii ));       
    end
    toc
    tic
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
    toc
    tic
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
    toc
    tic
    %% EM-KF
    epsilo = 0.00001;
    xe_EM = xe_ini;
    Pkk_EM = eye(n)*1;  
    

    
%     for ir = 1:size(r,1)
%         [pai(ir,:),R1(ir,:),R2(ir,:)] = EM_GMM2(r(ir,:));
%     end
%     R_pai = mean(pai);
%     RR1 = diag(R1);
%     RR2 = diag(R2);
    
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
%                         if ii > 200
%                 for  i = 1 : m 
%                     sigma = QQ1(i,i);
%                     Vik(i,1) = Q_pai(1) * exp(-1 * (eex(i) - 0)^2/(2*sigma)) / sqrt(2*pi*sigma); 
%                     sigma = QQ2(i,i);
%                     Vik(i,2) = Q_pai(2) * exp(-1 * (eex(i) - 0)^2/(2*sigma)) / sqrt(2*pi*sigma); 
%                     tt1(i,ii) = Vik(i,1) / (Vik(i,1) + Vik(i,2));
%                     tt2(i,ii) = Vik(i,2) / (Vik(i,1) + Vik(i,2));                   
%                     Cx(i,i) = ((tt1(i,ii) * (QQ2(i,i)) + tt2(i,ii) * (QQ1(i,i)))); 
% %                     Cx(i,i) = ((tt1(i,ii) * sqrt(QQ2(i,i)) + tt2(i,ii) * sqrt(QQ1(i,i)))/sqrt(Q(i,i))); 
%                 end    
%                         end
%             if ii == 500
%                 for ir = 1:size(r,1)
%                     [alphaK(ir,:), miuK(ir,:), sigmaK(ir,:)] = F_EM(ee_test(ir,end/2:end),2);
%                 end    
%                 R_pai = mean(alphaK);
%                 RR1 = diag(sigmaK(:,1));
%                 RR2 = diag(sigmaK(:,2)); 
% %                 Factor = R_pai(1) * (RR2) + R_pai(2) * (RR1);
%                 Miu1 = diag(miuK(:,1));
%                 Miu2 = diag(miuK(:,2));
%             end
%                 R_hat = Br * inv(Cy) * Br'
            PPKE = H*Pke;
            if ii > 50
                for  i = 1 : m 
                    sigma1 = RR1(i,i) + Pke(1,1);
                    Vik(i,1) = R_pai(1) * exp(-1 * (eey(i) - 0)^2/(2*sigma1)) / sqrt(2*pi*sigma1); 
                    sigma2 = RR2(i,i)+ Pke(1,1);
                    Vik(i,2) = R_pai(2) * exp(-1 * (eey(i) - 0)^2/(2*sigma2)) / sqrt(2*pi*sigma2); 
                    
                    t1(i,ii) = Vik(i,1) / (Vik(i,1) + Vik(i,2)); % 响应度
                    t2(i,ii) = Vik(i,2) / (Vik(i,1) + Vik(i,2));
                    
%                     tmp1 = t1(i,ii) * (RR1(i,i));
%                     tmp2 = t2(i,ii) * (RR2(i,i));
%                     Cy(i,i) = (tmp1 + tmp2);
                    Cy(i,i) = ((t1(i,ii) * (sigma2) + t2(i,ii) * (sigma1))) / R(i,i); 
%                     Cy(i,i) = ((t1(i,ii) * (sigma2) + t2(i,ii) * (sigma1))) ;
%                     Cy(i,i) = ((t1(i,ii) * sqrt(RR2(i,i)) + t2(i,ii) * sqrt(RR1(i,i)))/sqrt(R(i,i))); 
                   
%                     if Vik(i,1) == 0 && Vik(i,2) == 0
%                         Cy(i,i) = 1;
%                     end
                    
                    %% 
%                     Cy(i,i) = t1(ii) * (RR2(i,i)) + t2(ii) * (RR1(i,i));
                end    
%                 Cy = eye(m);
%                 tt1 = median(t1(:,ii));
%                 tt2 = median(t2(:,ii));
%                 for i = 1 : m
%                     tmp1 = tt1 * sqrt(RR2(i,i));
%                     tmp2 = tt2 * sqrt(RR1(i,i));
%                     Cy(i,i) = (tmp1 + tmp2) / sqrt(R(i,i));
%                 end
%                 vik1 = R_pai(1) * 1 / (2 * pi)^(m/2) / (abs(det(RR1)))^0.5 * exp(-0.5 * eey'*inv(RR1)*eey);
%                 vik2 = R_pai(2) * 1 / (2 * pi)^(m/2) / (abs(det(RR2)))^0.5 * exp(-0.5 * eey'*inv(RR2)*eey);
%                 t1(ii) = vik1 / (vik1 + vik2);
%                 t2(ii) = vik2 / (vik1 + vik2);
%                 Cy = ((t1(ii) * sqrt(RR2) + t2(ii) * sqrt(RR1)) ./ sqrt(R)); 
%                 Cy = ((t1(ii) * (RR2) + t2(ii) * (RR1)) ./ (R)); 
%                 Cy = diag(diag(Cy));
            end
            
            R_hat = Br * inv(Cy) * Br';
            
            R_hat = diag(diag(R_hat));
%             R_hat = Cy;
            Pke_hat = Bp * inv(Cx) * Bp';
            Gk_EM = Pke_hat * H' *inv(H*Pke_hat*H' + R_hat);
            xkk = xke + Gk_EM * (yy1 - H *  xke);
            xe_EM(:,ii) = xkk;
%             if (norm(xkk - past)/norm(past)) < epsilo
%                 xe_EM(:,ii) = xkk;
%                 break;
%             else
%                 past = xkk;          
%             end
%         end
        tmp1 = (eye(n) - Gk_EM*H);
        Pkk_EM = tmp1*Pke*tmp1' + Gk_EM*R*Gk_EM';
        Err_EM_KF(mm,ii) = (xkk - xx(:,ii))'*(xkk - xx(:,ii));
    end
    toc
    disp(mm)
end

figure, hold on;
plot(10*log10(mean(Err_MSE_KF)),'r');
% plot(10*log10(mean(Err_MCC_KF)),'g');
plot(10*log10(mean(Err_MEE_KF)),'b');
plot(10*log10(mean(Err_EM_KF)),'k');
xlabel('Iterations');ylabel('MSD(dB)');
