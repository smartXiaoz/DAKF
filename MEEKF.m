function [Pkk_MEEKF,X_MEEKF]=MEEKF(F,T,H,Q,R,Z,xe_MEEKF,Pkke_MEEKF,p,sigma)
theta=sigma;
t=1;
epsilo=0.1;
Q=T*Q*T';
xke = F * xe_MEEKF;
Pke = F * Pkke_MEEKF * F' + Q;
Bpk = chol(Pke)';
Brk = chol(R)';
Bkk = blkdiag(Bpk,Brk);
dik = inv(Bkk)*[xke;Z];
wik = inv(Bkk)*[eye(p);H];
xkk = xke;
while(t<1000)
    temp=xkk;
    t=t+1;
    ee = dik - wik * xkk;
    for i=1:length(ee)
        for j=1:length(ee)
            e=ee(i)-ee(j);
            phi_k(i,j)=exp(-e*e/(2*theta*theta));
        end
    end
    psi_k=diag(sum(phi_k));
    Omiga=phi_k^2+psi_k^2;
    Px=(inv(Bpk))'*Omiga(1:p,1:p)*inv(Bpk);
    Pxy=(inv(Brk))'*Omiga(p+1:end,1:p)*inv(Bpk);
    Pyx=(inv(Bpk))'*Omiga(1:p,p+1:end)*inv(Brk);
    Py=(inv(Brk))'*Omiga(p+1:end,p+1:end)*inv(Brk);
    Gk_MEEKF = inv(wik'*Omiga*wik)*(Pyx+H'*Py);
    xkk = xke + Gk_MEEKF*(Z-H * xke);
    xe_MEEKF = xkk;
    
    if(abs(xe_MEEKF-temp)/abs(temp)<=epsilo)
        break;
    else
        t=t-1;
    end
end
X_MEEKF=xe_MEEKF;
Pkke_MEEKF = (eye(p) - Gk_MEEKF * H) * Pke * (eye(p) - Gk_MEEKF * H)'+ Gk_MEEKF * R * Gk_MEEKF';
Pkk_MEEKF=Pkke_MEEKF;
end
