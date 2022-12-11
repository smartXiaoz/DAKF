function [xeii_t,Pk_MEE,e,Phi_k,Psi_k,lambda,Ke] = MEE_KF_F(F,xx1,Pk_MEE,H,yy1,Q,R,sigma_MEE)

m = length(yy1);
n = length(xx1);
epsilon = 1e-6;
Pke = F * Pk_MEE * F' + Q;
xeii_t0 = F * xx1; 
xeii_t = xeii_t0 ; 
xeii_t1 = xeii_t0 ;
%     theta_all = [Pke,zeros(n,m);zeros(m,n),R];
%     theta  = (chol(theta_all))';
%     theta_p = theta(1:n,1:n);     theta_r = theta((n+1):(n+m),(n+1):(n+m));
theta_p =(chol(Pke))' ;     theta_r =(chol(R))' ;
theta = [theta_p zeros(n, m); zeros(m, n)  theta_r];
D  = inv(theta) * [xeii_t0 ; yy1];  W  = inv(theta) * [eye(n); H];
Yita = yy1 - H* xeii_t0;
A_Pke_R = H * Pke * H' + R;



qq = m+n;
for b = 1 : 20
    xeii_t1 = xeii_t;
    %error
    e = D - W * xeii_t1;
    
    Phi_k = zeros(n+m,n+m);
    Psi_k = zeros(n+m,n+m);
    for kk = 1:qq
        for jj = 1:qq
            eij = e(jj) - e(kk);
            Phi_k(kk,jj) = exp(-eij^2/(2*sigma_MEE^2)) ;
        end
    end
    Psi_k = diag( sum( Phi_k ) );
    lambda =Psi_k^2 + Phi_k^2;
    %               lambda =1*Psi_k - 1.0*Phi_k;
    %              lambda =Psi_k - diag(diag(Phi_k));
    lambda_x = lambda(1:n,1:n);
    lambda_y = lambda(n+1:n+m,n+1:n+m);
    lambda_xy = lambda(n+1:n+m,1:n);
    lambda_yx = lambda(1:n,n+1:n+m);
    %         lambda_xy =0; lambda_yx =0  ;
    
    %% ¶¨µãµü´ú
    Pxx =  inv(theta_p)' * lambda_x * inv(theta_p) ;
    Pxy = inv(theta_r)' * lambda_xy * inv(theta_p) ;
    Pyx =  inv(theta_p)' * lambda_yx * inv(theta_r) ;
    Pyy =  inv(theta_r)' * lambda_y * inv(theta_r) ;
    tmp = inv(W'*lambda*W);
    Ke = tmp* (Pyx + H' * Pyy);
    
    %% posterior mean update
    xeii_t = xeii_t0 + Ke * (yy1 - H * xeii_t0);
    %% loop counter
    
    if ( norm( xeii_t - xeii_t1 )<= epsilon )
        break;
    end
end

%% optimal estimated state of current time step

%% estimated state covariance of current time step
Pk_MEE = ( eye(n) - Ke * H ) * Pke * ( eye(n)-Ke * H )'+Ke * R * Ke';
