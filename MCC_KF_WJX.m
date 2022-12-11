function [xx3,Pk_MCC,b,C,Ke_MCC,M_MCC] = MCC_KF_F(F,xx3,Pk_MCC,H,yy3,Q,R,sigma_MCC)
%% mckf maximum correntropy Kalman filter
% [x,P,b]=mckf(F,x,P,H,z,Q,R)
% Inputs:  F: state transition matrix
%          x: optimal estimated state of last time step
%          P: estimated state covariance of last time step
%          H: observation matrix
%          z: measurement of current time step
%          Q: covariance matrix of process noise
%          R: covariance matrix of measurement noise
%          b: iteration number
% Outputs: x: optimal estimated state of current time step
%          P: estimated state covariance of current time step

%% 定点迭代参数配置
n  = numel(xx3);   %number of states
m  = numel(yy3);   %number of measurements
epsilon = 1e-6;
xeii_t0 = F * xx3;
xxii_t = xeii_t0+100 ;
% xeii_tt =xeii_t0 +  ones(1);%
xeii_tt =xeii_t0;
%xeii_tt = eye(1);
Pke = F * Pk_MCC * F' + Q;   %priori estimated state covariance
% B_all = [Pke,zeros(n,m);zeros(m,n),R];
% B  = (chol(B_all))';
theta_p =(chol(Pke))' ;     
theta_r =(chol(R))' ;
theta = inv([theta_p zeros(n, m); zeros(m, n)  theta_r]);
D  = theta* [xeii_t0 ; yy3];
W  = theta* [eye(n); H];
b = 0;
theta_p =theta(1:n,1:n) ;     
theta_r =theta((n+1):(n+m),(n+1):(n+m)) ;
%% Now iterate
while (norm((xeii_tt - xxii_t)/norm(xxii_t))>=epsilon && (b<=10000))
    xxii_t = xeii_tt;
    %error
    e = D - W * xxii_t;
    % Correntropy calculation
    for a = 1:n+m
        G(a) = exp(-(e(a)^2)/(2 * sigma_MCC^2));
%         G(a) = -e(a)^(2);
        %G(a) = 1;
    end  
    % Diagonalise Correntropies
    C   = diag(G);  
    
    % 'Measured' covariance update
    Pee = theta_p' *(C(1:n,1:n)) * theta_p;
    Re_all=theta_r'*(C((n+1):(n+m),(n+1):(n+m)))* theta_r;
    % Propagated covariance update
    Re = H'*Re_all*H;  
    % Kalman gain
    M_MCC = inv(Pee+Re);
    % posterior mean update
    xeii_tt = M_MCC*(Pee*xeii_t0+ H'*Re_all*yy3); 

    % loop counter
    b = b + 1;
end
Ke_MCC = M_MCC * H' * Re_all;
%% optimal estimated state of current time step
xx3 = xeii_tt;

%% estimated state covariance of current time step
% Pk_MCC = M_MCC;
Pk_MCC = ( eye(n) - Ke_MCC * H ) * Pke * ( eye(n)-Ke_MCC * H )' ...
+  Ke_MCC * R * Ke_MCC';