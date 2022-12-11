function Err_MEE_TH = fun_steady_state_error_KF(H,F,Q,Ke1,KeR_MEE);
%% x: n*1; y: m*1; Q: n*n; Ke1: n*m; r_est2: m*m;

n = size(Q,1);
Yk = ( eye(n) - Ke1 * H ) * F;
Zk = ( eye(n) - Ke1 * H ) * Q * ( eye(n)-Ke1 * H )'+KeR_MEE;
Zk = Zk(:);
Vec_e = inv( eye(n^2) - kron(Yk,Yk) ) * Zk;
Mat_e = zeros(n);
Mat_e(:) = Vec_e(:);
Err_MEE_TH  = trace(Mat_e);