function [pai,R1,R2]=EM_GMM2(v,mu1,mu2)
q=size(v,1);
Len=size(v,2);
R=v*v'/Len;
R_EM(:,:,1)=R;
R_EM(:,:,2)=R/2;
pai(1)=0.5;
pai(2)=0.5;
mu=zeros(q,2);
mu=[mu1,mu2;mu1,mu2];
dabul=zeros(2,Len);
vpa(dabul,16);%将数据精度调整至小数点后50位
EM=1;
while(EM<30)
    temp=mu;
    %Expectation
    for i=1:Len
        for tar=1:2
            dzz=v(:,i)-mu(:,tar);         
%             ee_2 = dik_2 - wik_2 * xkk_2;
%             Gee_2 = 0.9*exp(-ee_2.^2/2/sigma^2)+0.1*tanh(0.2*dik_2.*(wik_2 * xkk_2)+5);%IMCKF核
%             dabul(tar,i)= 0.9*pai(tar)*1/(2*pi)^(length(mu(:,tar))/2)/(abs(det(R_EM(:,:,tar))))^0.5*exp(-0.5*dzz'*inv(R_EM(:,:,tar))*dzz)+...
%                 0.1*tanh(0.2*v(:,i)'*mu(:,tar)+5);

            dabul(tar,i)= pai(tar)*1/(2*pi)^(length(mu(:,tar))/2)/(abs(det(R_EM(:,:,tar))))^0.5*exp(-0.5*dzz'*inv(R_EM(:,:,tar))*dzz);
        end
        nfac=sum(dabul(:,i));
        if(isnan(nfac)==1)
%             dabul(:,i)= 1/2;
            dabul(:,i)=pai;
        else
            dabul(:,i)= dabul(:,i)/nfac;  
        end

    end
        %Max
        for tar=1:2
                pai(tar)=sum(dabul(tar,:))/Len;
                sum_dabul_meas=zeros(q,1);
                sum_dabul_R=zeros(q,q);
                for mea =1:Len
                    sum_dabul_meas=sum_dabul_meas+dabul(tar,mea)*v(:,mea);
                    sum_dabul_R=sum_dabul_R+dabul(tar,mea)*(v(:,mea)-mu(:,tar))*(v(:,mea)-mu(:,tar))';
                end
                mu(:,tar)=sum_dabul_meas/sum(dabul(tar,:));
                R_EM(:,:,tar)=sum_dabul_R/sum(dabul(tar,:));
        end
%         pai_EM(1,EM)=pai(1);
%         pai_EM(2,EM)=pai(2);
        if(sum(sum(abs(mu-temp)))<0.0001)
            break;
        else
        EM=EM+1;
        end
end
R1=R_EM(:,:,1);
R2=R_EM(:,:,2);
end