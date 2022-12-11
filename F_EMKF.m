function [Pkk_MCKF3,X_MCKF3]=F_EMKF(F,T,H,Q,R,Z,xe_MCKF,Pkk_MCKF,p,RR1,RR2,R_pai)
    sigma=5;
    t=1;
    epsilo=0.1;
    Q=T*Q*T';
    xke = F * xe_MCKF;
    Pke = F * Pkk_MCKF * F' + Q;
    Bpk = chol(Pke)';
    Brk = chol(R)';
    Bkk = blkdiag(Bpk,Brk);
    dik = inv(Bkk)*[xke;Z];
    wik = [inv(Bpk);inv(Brk)* H];
    wik = inv(Bkk)*[eye(p);H];
    n = size(Bpk,1);
    m = size(Brk,1);
    xkk = xke;
%   while(t<1000)
      temp=xkk;
       t=t+1;
        ee = dik - wik * xkk;
        Gee = exp(-ee.^2/2/sigma^2);
        Cx = eye(p);
%         Cy = diag(Gee(p+1:end));
%         Pke_hat_inv = inv(Bpk)' * (Cx) * inv(Bpk);
%         R_hat_inv = inv(Brk)' * (Cy) * inv(Brk);
%         Gk_DMCKF = inv(H'*R_hat_inv*H+Pke_hat_inv)*H'*R_hat_inv;%DMCKF
        eey = Z - H * xke;
        PP1 = H*Pke*H' + RR1;
        PP2 = H*Pke*H' + RR2;
        vik1 = R_pai(1) * 1 / (2 * pi)^(m/2) / (abs(det(PP1)))^0.5 * exp(-0.5 * eey'*inv(PP1)*eey);
        vik2 = R_pai(2) * 1 / (2 * pi)^(m/2) / (abs(det(PP2)))^0.5 * exp(-0.5 * eey'*inv(PP2)*eey);
        sumaa = (vik1 + vik2);
        if sumaa == 0
            t1 = 1;
            t2 = 0;
        else
            t1 = vik1 / sumaa;
            t2 = vik2 / sumaa;
        end
%        if abs(sum(eey)) > 3 * sum(diag(PP2.^0.5))
%             t1 = 1;
%             t2 = 0;
%         else
%              t1 = 0;
%              t2 = 1;
%         end
                
        Cy = ((t1 * RR2.^0.5 + t2 * RR1.^0.5) ./ R.^0.5); 
%         for i = 1 : m
%             if abs(eey(i)) > 3 * sqrt(PP2(i,i))
%                 Cy(i,i) = sqrt(RR2(i,i))/ sqrt(R(i,i));              
%             else
%                 Cy(i,i) = sqrt(RR1(i,i))/ sqrt(R(i,i)); 
%             end
%         end
        
%         Cy = (t1 * RR2 + t2 * RR1) ./ R; 
        Cy = diag(diag(Cy));
        Pke_hat = Bpk * inv(Cx) * Bpk';
        R_hat = Brk * inv(Cy) * Brk';
        Gk_DMCKF = Pke_hat * H' * inv(H * Pke_hat * H'+R_hat);%MCKF
        
        xkk = xke + Gk_DMCKF*(Z-H * xke);
        
       xe_MCKF = xkk;
    
%       if(abs(xe_MCKF-temp)/abs(temp)<=epsilo)
%           break;
%       else
%       t=t-1;
%       end
%   end
    X_MCKF3=xe_MCKF;
    Pkk_MCKF = (eye(p) - Gk_DMCKF * H) * Pke * (eye(p) - Gk_DMCKF * H)'+ Gk_DMCKF * R * Gk_DMCKF';
    Pkk_MCKF3=Pkk_MCKF;
    
% end
