function [Pkk_MCKF3,X_MCKF3]=MCKF3(F,T,H,Q,R,Z,xe_MCKF,Pkk_MCKF,p,sigma)
    
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
    
    xkk = xke;
  while(t<10)
      temp=xkk;
       t=t+1;
        ee = dik - wik * xkk;
        Gee = exp(-ee.^2/2/sigma^2);
        Cx = diag(Gee(1:p));
        Cy = diag(Gee(p+1:end));
%         Pke_hat_inv = inv(Bpk)' * (Cx) * inv(Bpk);
%         R_hat_inv = inv(Brk)' * (Cy) * inv(Brk);
%         Gk_DMCKF = inv(H'*R_hat_inv*H+Pke_hat_inv)*H'*R_hat_inv;%DMCKF
%       
        Pke_hat = Bpk * inv(Cx) * Bpk';
        R_hat = Brk * inv(Cy) * Brk';
        Gk_DMCKF = Pke_hat * H' * inv(H * Pke_hat * H'+R_hat);%MCKF
        
        xkk = xke + Gk_DMCKF*(Z-H * xke);
        
       xe_MCKF = xkk;
    
      if(abs(xe_MCKF-temp)/abs(temp)<=epsilo)
          break;
      else
      t=t-1;
      end
  end
    X_MCKF3=xe_MCKF;
    Pkk_MCKF = (eye(p) - Gk_DMCKF * H) * Pke * (eye(p) - Gk_DMCKF * H)'+ Gk_DMCKF * R * Gk_DMCKF';
    Pkk_MCKF3=Pkk_MCKF;
    
end
