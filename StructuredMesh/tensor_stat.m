close all; 
clear; 
clc
double precision;
format long;

%% Load data
load 'data_stats.mat';
load 'angles.mat';
load 'geo.mat';
pipe=struct;

%%re-order x,y
% npoints=length(x);
% k=0;
% for i=1:nr
%     for j=1:nth
%         k=k+1;
%         X(i,j)=x(k);
%         Y(i,j)=y(k);
%     end
% end
% k=0;
% for j=1:nth
%     for i=1:nr
%         k=k+1;
%         x1(k,1)=X(i,j);
%         y1(k,1)=Y(i,j);
%     end
% end
% 
% for i = 1:length(zn)
%     f=length(x1);
%     x(1+(i-1)*f:i*f) = x1;
%     y(1+(i-1)*f:i*f) = y1;
% end
% % 
% % transpose the imported data tensors
% U=U'; V=V'; W=W';
% uu=uu'; uv=uv'; uw=uw'; vv=vv'; vw=vw'; ww=ww';
% dUdx=dUdx'; dUdy=dUdy'; dUdz=dUdz';
% dVdx=dVdx'; dVdy=dVdy'; dVdz=dVdz';
% dWdx=dWdx'; dWdy=dWdy'; dWdz=dWdz';
% Pxx=Pxx'; Pxy=Pxy'; Pxz=Pxz'; Pyy=Pyy'; Pyz=Pyz'; Pzz=Pzz'; 
% Dxx=Dxx'; Dxy=Dxy'; Dxz=Dxz'; Dyy=Dyy'; Dyz=Dyz'; Dzz=Dzz'; 
% Cxx=Cxx'; Cxy=Cxy'; Cxz=Cxz'; Cyy=Cyy'; Cyz=Cyz'; Czz=Czz'; 
% Txx=Txx'; Txy=Txy'; Txz=Txz'; Tyy=Tyy'; Tyz=Tyz'; Tzz=Tzz'; 
% VDxx=VDxx'; VDxy=VDxy'; VDxz=VDxz'; VDyy=VDyy'; VDyz=VDyz'; VDzz=VDzz'; 
% PTxx=PTxx'; PTxy=PTxy'; PTxz=PTxz'; PTyy=PTyy'; PTyz=PTyz'; PTzz=PTzz'; 
% PSxx=PSxx'; PSxy=PSxy'; PSxz=PSxz'; PSyy=PSyy'; PSyz=PSyz'; PSzz=PSzz'; 
% uuu=uuu'; uvv=uvv'; uuw=uuw';
% uuv=uuv';           uvw=uvw';
% uuw=uuw';           uww=uww';
%           vvv=vvv'; vvw=vvw';
%                     vww=vww';
%                     www=www'; 
% % duuudx = duuudx';, 'duuvdx','duuwdx', 'duvvdx', 'duvwdx', 'duwwdx', ...
% %       'duuvdy', 'duvvdy', 'duvwdy', 'dvvvdy', 'dvvwdy', 'dvwwdy', ...
% %       'duuwdz', 'duvwdz', 'duwwdz', 'dvvwdz', 'dvwwdz', 'dwwwdz');                    
%% Rotation to cylindrical coord.

rho=1;
n2=0;
n3=0;
    
for k=1:nz
    for i=1:nth
        
        n1=n2+1;
        n2=nr*i+n3;
        
        pipe.th1(i,:,k) = angle(n1:n2);
        pipe.r1(i,:,k)  = sqrt(x(n1:n2).^2+y(n1:n2).^2);
        pipe.z1(i,:,k)  = z(n1:n2);
        
        R = [cos(angle(i)) sin(angle(i)) 0;
            -sin(angle(i)) cos(angle(i)) 0;
                 0             0            1];
             
        for jj=n1:n2
            j=jj-n1+1;
            
            %Mean velocities. Tensors of Rank 1.
            prod = R*[U(jj); V(jj); W(jj)];       
            
            pipe.Ur(i,j,k) = prod(1);
            pipe.Ut(i,j,k) = prod(2);
            pipe.Uz(i,j,k) = prod(3); 
        
            %Reynolds stress tensor. Tensors of Rank 2.
            S = [uu(jj) uv(jj) uw(jj); 
                 uv(jj) vv(jj) vw(jj); 
                 uw(jj) vw(jj) ww(jj)];
     
            prod = R*S*R';

            pipe.urur(i,j,k) = prod(1,1);
            pipe.urut(i,j,k) = prod(1,2);
            pipe.uruz(i,j,k) = prod(1,3);
            pipe.utut(i,j,k) = prod(2,2);
            pipe.uzuz(i,j,k) = prod(3,3);
            pipe.utuz(i,j,k) = prod(2,3);
            
            %Velocity gradient tensor. Tensor of Rank 2.     
            prod = R*[dUdx(jj) dUdy(jj) dUdz(jj);
                     dVdx(jj) dVdy(jj) dVdz(jj);
                     dWdx(jj) dWdy(jj) dWdz(jj)]*R';
        
            pipe.dUrdr(i,j,k)  = prod(1,1);
            pipe.dUrdt(i,j,k)  = prod(1,2);
            pipe.dUrdz(i,j,k)  = prod(1,3);
            pipe.dUtdr(i,j,k)  = prod(2,1);
            pipe.dUtdt(i,j,k)  = prod(2,2);
            pipe.dUtdz(i,j,k)  = prod(2,3);
            pipe.dUzdr(i,j,k)  = prod(3,1);
            pipe.dUzdt(i,j,k)  = prod(3,2);   
            pipe.dUzdz(i,j,k)  = prod(3,3);
 
            %Production tensor. Tensor of Rank 2.
            S = [Pxx(jj) Pxy(jj) Pxz(jj); 
                 Pxy(jj) Pyy(jj) Pyz(jj); 
                 Pxz(jj) Pyz(jj) Pzz(jj)];

            prod=R*S*R';        

            pipe.Prr(i,j,k) = prod(1,1);
            pipe.Ptt(i,j,k) = prod(2,2);
            pipe.Pzz(i,j,k) = prod(3,3);
            pipe.Prt(i,j,k) = prod(1,2);
            pipe.Prz(i,j,k) = prod(1,3);
            pipe.Ptz(i,j,k) = prod(2,3);
 
            %Dissipation tensor. Tensor of Rank 2.
            S = [Dxx(jj) Dxy(jj) Dxz(jj); 
                 Dxy(jj) Dyy(jj) Dyz(jj);  
                 Dxz(jj) Dyz(jj) Dzz(jj)];

            prod=R*S*R';

            pipe.Drr(i,j,k) = prod(1,1);
            pipe.Dtt(i,j,k) = prod(2,2);
            pipe.Dzz(i,j,k) = prod(3,3);
            pipe.Drt(i,j,k) = prod(1,2);
            pipe.Drz(i,j,k) = prod(1,3);
            pipe.Dtz(i,j,k) = prod(2,3);
 
            %Mean convection tensor. Tensor of Rank 2.
            S = [Cxx(jj) Cxy(jj) Cxz(jj); 
                 Cxy(jj) Cyy(jj) Cyz(jj);  
                 Cxz(jj) Cyz(jj) Czz(jj)];

            prod=R*S*R';   

            pipe.Crr(i,j,k) = prod(1,1);
            pipe.Ctt(i,j,k) = prod(2,2);
            pipe.Czz(i,j,k) = prod(3,3);
            pipe.Crt(i,j,k) = prod(1,2);
            pipe.Crz(i,j,k) = prod(1,3);
            pipe.Ctz(i,j,k) = prod(2,3);
    
            %Turbulent transport tensor. Tensor of Rank 2.
            S = [Txx(jj) Txy(jj) Txz(jj); 
                 Txy(jj) Tyy(jj) Tyz(jj);  
                 Txz(jj) Tyz(jj) Tzz(jj)];

            prod=R*S*R'; 

            pipe.Trr(i,j,k) = prod(1,1);
            pipe.Ttt(i,j,k) = prod(2,2);
            pipe.Tzz(i,j,k) = prod(3,3);
            pipe.Trt(i,j,k) = prod(1,2);
            pipe.Trz(i,j,k) = prod(1,3);
            pipe.Ttz(i,j,k) = prod(2,3);

            %Viscous diffusion tensor. Tensor of Rank 2.
            S = [VDxx(jj) VDxy(jj) VDxz(jj); 
                 VDxy(jj) VDyy(jj) VDyz(jj);  
                 VDxz(jj) VDyz(jj) VDzz(jj)];

            prod=R*S*R';    

            pipe.VDrr(i,j,k) = prod(1,1);
            pipe.VDtt(i,j,k) = prod(2,2);
            pipe.VDzz(i,j,k) = prod(3,3);
            pipe.VDrt(i,j,k) = prod(1,2);
            pipe.VDrz(i,j,k) = prod(1,3);
            pipe.VDtz(i,j,k) = prod(2,3);

            %Pressure transport tensor. Tensor of Rank 2.   
            S = [PTxx(jj) PTxy(jj) PTxz(jj); 
                 PTxy(jj) PTyy(jj) PTyz(jj);  
                 PTxz(jj) PTyz(jj) PTzz(jj)];

            prod=R*S*R';    

%             pipe.PTrr(i,j,k) = -2/rho*prod(1,1);
%             pipe.PTtt(i,j,k) = -2/rho*prod(2,2);
%             pipe.PTzz(i,j,k) = -2/rho*prod(3,3);
%             pipe.PTrt(i,j,k) = -1/rho*prod(1,2);
%             pipe.PTrz(i,j,k) = -1/rho*prod(1,3);
%             pipe.PTtz(i,j,k) = -1/rho*prod(2,3); 

            pipe.PTrr(i,j,k) = prod(1,1);
            pipe.PTtt(i,j,k) = prod(2,2);
            pipe.PTzz(i,j,k) = prod(3,3);
            pipe.PTrt(i,j,k) = prod(1,2);
            pipe.PTrz(i,j,k) = prod(1,3);
            pipe.PTtz(i,j,k) = prod(2,3);
            
            %Pressure strain tensor. Tensor of Rank 2.
            S = [PSxx(jj) PSxy(jj) PSxz(jj); 
                 PSxy(jj) PSyy(jj) PSyz(jj);  
                 PSxz(jj) PSyz(jj) PSzz(jj)];

            prod=R*S*R';    

%             pipe.PSrr(i,j,k) = -2/rho*prod(1,1);
%             pipe.PStt(i,j,k) = -2/rho*prod(2,2);
%             pipe.PSzz(i,j,k) = -2/rho*prod(3,3);
%             pipe.PSrt(i,j,k) = -1/rho*prod(1,2);
%             pipe.PSrz(i,j,k) = -1/rho*prod(1,3);
%             pipe.PStz(i,j,k) = -1/rho*prod(2,3); 

            pipe.PSrr(i,j,k) = prod(1,1);
            pipe.PStt(i,j,k) = prod(2,2);
            pipe.PSzz(i,j,k) = prod(3,3);
            pipe.PSrt(i,j,k) = prod(1,2);
            pipe.PSrz(i,j,k) = prod(1,3);
            pipe.PStz(i,j,k) = prod(2,3);
    
           %Budget for each component of the Reynolds stress tensor 
%            %Without mean convection
%             pipe.Srr(i,j,k) = pipe.Prr(i,j,k)+pipe.Drr(i,j,k)+pipe.Trr(i,j,k)+pipe.VDrr(i,j,k)+pipe.Pirr(i,j,k); 
%             pipe.Stt(i,j,k) = pipe.Ptt(i,j,k)+pipe.Dtt(i,j,k)+pipe.Ttt(i,j,k)+pipe.VDtt(i,j,k)+pipe.Pitt(i,j,k);
%             pipe.Szz(i,j,k) = pipe.Pzz(i,j,k)+pipe.Dzz(i,j,k)+pipe.Tzz(i,j,k)+pipe.VDzz(i,j,k)+pipe.Pizz(i,j,k);
%             pipe.Srt(i,j,k) = pipe.Prt(i,j,k)+pipe.Drt(i,j,k)+pipe.Trt(i,j,k)+pipe.VDrt(i,j,k)+pipe.Pirt(i,j,k); 
%             pipe.Srz(i,j,k) = pipe.Prz(i,j,k)+pipe.Drz(i,j,k)+pipe.Trz(i,j,k)+pipe.VDrz(i,j,k)+pipe.Pirz(i,j,k);
%             pipe.Stz(i,j,k) = pipe.Ptz(i,j,k)+pipe.Dtz(i,j,k)+pipe.Ttz(i,j,k)+pipe.VDtz(i,j,k)+pipe.Pitz(i,j,k);
%            %With mean convection
%             pipe.Scrr(i,j,k) = pipe.Prr(i,j,k)+pipe.Drr(i,j,k)+pipe.Trr(i,j,k)+pipe.VDrr(i,j,k)+pipe.Pirr(i,j,k)-pipe.Crr(i,j,k); 
%             pipe.Sctt(i,j,k) = pipe.Ptt(i,j,k)+pipe.Dtt(i,j,k)+pipe.Ttt(i,j,k)+pipe.VDtt(i,j,k)+pipe.Pitt(i,j,k)-pipe.Ctt(i,j,k);
%             pipe.Sczz(i,j,k) = pipe.Pzz(i,j,k)+pipe.Dzz(i,j,k)+pipe.Tzz(i,j,k)+pipe.VDzz(i,j,k)+pipe.Pizz(i,j,k)-pipe.Czz(i,j,k);
%             pipe.Scrt(i,j,k) = pipe.Prt(i,j,k)+pipe.Drt(i,j,k)+pipe.Trt(i,j,k)+pipe.VDrt(i,j,k)+pipe.Pirt(i,j,k)-pipe.Crt(i,j,k); 
%             pipe.Scrz(i,j,k) = pipe.Prz(i,j,k)+pipe.Drz(i,j,k)+pipe.Trz(i,j,k)+pipe.VDrz(i,j,k)+pipe.Pirz(i,j,k)-pipe.Crz(i,j,k);
%             pipe.Sctz(i,j,k) = pipe.Ptz(i,j,k)+pipe.Dtz(i,j,k)+pipe.Ttz(i,j,k)+pipe.VDtz(i,j,k)+pipe.Pitz(i,j,k)-pipe.Ctz(i,j,k);

            %Skewness tensor. Tensor of Rank 3.    
            R3_tensor(:,:,1) = [uuu(jj) uvv(jj) uuw(jj);...
                                uuv(jj) uvv(jj) uvw(jj);...
                                uuw(jj) uvw(jj) uww(jj)];
            R3_tensor(:,:,2) = [uuv(jj) uvv(jj) uvw(jj);...
                                uvv(jj) vvv(jj) vvw(jj);...
                                uvw(jj) vvw(jj) vww(jj)];
            R3_tensor(:,:,3) = [uuw(jj) uvw(jj) uww(jj);...
                                uvw(jj) vvw(jj) vww(jj);...
                                uww(jj) vww(jj) www(jj)];

            aabc(1:3,1:3,1:3)=0;
            adef=R3_tensor(:,:,:);
            for aa=1:3
            for bb=1:3    
            for cc=1:3    
            for dd=1:3
            for ee=1:3    
            for ff=1:3
                aabc(aa,bb,cc)=aabc(aa,bb,cc)+R(aa,dd)*R(bb,ee) ...
                    *R(cc,ff)*adef(dd,ee,ff);
            end    
            end
            end
            end    
            end
            end

            pipe.ururur(i,j,k) = aabc(1,1,1);
            pipe.ututut(i,j,k) = aabc(2,2,2);
            pipe.uzuzuz(i,j,k) = aabc(3,3,3);
            pipe.ururut(i,j,k) = aabc(1,2,1);
            pipe.ururuz(i,j,k) = aabc(1,3,1);
            pipe.urutut(i,j,k) = aabc(2,2,1);
            pipe.ututuz(i,j,k) = aabc(2,3,2);
            pipe.uruzuz(i,j,k) = aabc(3,3,1);
            pipe.utuzuz(i,j,k) = aabc(3,3,2);
            pipe.urutuz(i,j,k) = aabc(2,3,1);    
          
        end       
    end    
    n3=n2;
end

%% Azimuthal Average

%clearvars -except pipe nz nr nu yplus zn x y
pipe2=struct;

for k=1:nz
    
    for j=1:nr
        
        %Mean velocities. Tensors of Rank 1.        
        
        pipe2.Ur(j,k) = mean(pipe.Ur(:,j,k));
        pipe2.Ut(j,k) = mean(pipe.Ut(:,j,k));       
        pipe2.Uz(j,k) = mean(pipe.Uz(:,j,k));

        %Reynolds stress tensor. Tensors of Rank 2.
        
        pipe2.urur(j,k) = mean(pipe.urur(:,j,k));
        pipe2.urut(j,k) = mean(pipe.urut(:,j,k));
        pipe2.uruz(j,k) = mean(pipe.uruz(:,j,k));
        pipe2.utut(j,k) = mean(pipe.utut(:,j,k));
        pipe2.uzuz(j,k) = mean(pipe.uzuz(:,j,k));
        pipe2.utuz(j,k) = mean(pipe.utuz(:,j,k));
        
        %Velocity gradient tensor. Tensor of Rank 2.
        
        pipe2.dUrdr(j,k) = mean(pipe.dUrdr(:,j,k));
        pipe2.dUrdt(j,k) = mean(pipe.dUrdt(:,j,k));
        pipe2.dUrdz(j,k) = mean(pipe.dUrdz(:,j,k));
        pipe2.dUtdr(j,k) = mean(pipe.dUtdr(:,j,k));
        pipe2.dUtdt(j,k) = mean(pipe.dUtdt(:,j,k));
        pipe2.dUtdz(j,k) = mean(pipe.dUtdz(:,j,k));
        pipe2.dUzdr(j,k) = mean(pipe.dUzdr(:,j,k));
        pipe2.dUzdt(j,k) = mean(pipe.dUzdt(:,j,k)); 
        pipe2.dUzdz(j,k) = mean(pipe.dUzdz(:,j,k));    
         
        %Production tensor. Tensor of Rank 2.     

        pipe2.Prr(j,k) = mean(pipe.Prr(:,j,k));
        pipe2.Ptt(j,k) = mean(pipe.Ptt(:,j,k));
        pipe2.Pzz(j,k) = mean(pipe.Pzz(:,j,k));
        pipe2.Prt(j,k) = mean(pipe.Prt(:,j,k));
        pipe2.Prz(j,k) = mean(pipe.Prz(:,j,k));
        pipe2.Ptz(j,k) = mean(pipe.Ptz(:,j,k));

        %Dissipation tensor. Tensor of Rank 2.
       
        pipe2.Drr(j,k) = mean(pipe.Drr(:,j,k));
        pipe2.Dtt(j,k) = mean(pipe.Dtt(:,j,k));
        pipe2.Dzz(j,k) = mean(pipe.Dzz(:,j,k));
        pipe2.Drt(j,k) = mean(pipe.Drt(:,j,k));
        pipe2.Drz(j,k) = mean(pipe.Drz(:,j,k));
        pipe2.Dtz(j,k) = mean(pipe.Dtz(:,j,k));
  
        %Mean convection tensor. Tensor of Rank 2.
        
        pipe2.Crr(j,k) = mean(pipe.Crr(:,j,k));
        pipe2.Ctt(j,k) = mean(pipe.Ctt(:,j,k));
        pipe2.Czz(j,k) = mean(pipe.Czz(:,j,k));
        pipe2.Crt(j,k) = mean(pipe.Crt(:,j,k));
        pipe2.Crz(j,k) = mean(pipe.Crz(:,j,k));
        pipe2.Ctz(j,k) = mean(pipe.Ctz(:,j,k));
  
        %Turbulent transport tensor. Tensor of Rank 2.

        pipe2.Trr(j,k) = mean(pipe.Trr(:,j,k));
        pipe2.Ttt(j,k) = mean(pipe.Ttt(:,j,k));
        pipe2.Tzz(j,k) = mean(pipe.Tzz(:,j,k));
        pipe2.Trt(j,k) = mean(pipe.Trt(:,j,k));
        pipe2.Trz(j,k) = mean(pipe.Trz(:,j,k));
        pipe2.Ttz(j,k) = mean(pipe.Ttz(:,j,k));

        %Viscous diffusion tensor. Tensor of Rank 2.
       
        pipe2.VDrr(j,k) = mean(pipe.VDrr(:,j,k));
        pipe2.VDtt(j,k) = mean(pipe.VDtt(:,j,k));
        pipe2.VDzz(j,k) = mean(pipe.VDzz(:,j,k));
        pipe2.VDrt(j,k) = mean(pipe.VDrt(:,j,k));
        pipe2.VDrz(j,k) = mean(pipe.VDrz(:,j,k));
        pipe2.VDtz(j,k) = mean(pipe.VDtz(:,j,k));
 
        %Pressure transport tensor. Tensor of Rank 2.   
        pipe2.PTrr(j,k) = mean(pipe.PTrr(:,j,k));
        pipe2.PTtt(j,k) = mean(pipe.PTtt(:,j,k));
        pipe2.PTzz(j,k) = mean(pipe.PTzz(:,j,k));
        pipe2.PTrt(j,k) = mean(pipe.PTrt(:,j,k));
        pipe2.PTrz(j,k) = mean(pipe.PTrz(:,j,k));
        pipe2.PTtz(j,k) = mean(pipe.PTtz(:,j,k));

        %Pressure strain tensor. Tensor of Rank 2.
        pipe2.PSrr(j,k) = mean(pipe.PSrr(:,j,k));
        pipe2.PStt(j,k) = mean(pipe.PStt(:,j,k));
        pipe2.PSzz(j,k) = mean(pipe.PSzz(:,j,k));
        pipe2.PSrt(j,k) = mean(pipe.PSrt(:,j,k));
        pipe2.PSrz(j,k) = mean(pipe.PSrz(:,j,k));
        pipe2.PStz(j,k) = mean(pipe.PStz(:,j,k));

%         
%         %Budget for each component of the Reynolds stress tensor 
%         %Without mean convection
%         pipe2.Srr(j,k) = mean(pipe.Srr(:,j,k));
%         pipe2.Stt(j,k) = mean(pipe.Stt(:,j,k));
%         pipe2.Szz(j,k) = mean(pipe.Szz(:,j,k));
%         pipe2.Srt(j,k) = mean(pipe.Srt(:,j,k));
%         pipe2.Srz(j,k) = mean(pipe.Srz(:,j,k));
%         pipe2.Stz(j,k) = mean(pipe.Stz(:,j,k));
%         %With mean convection
%         pipe2.Scrr(j,k) = mean(pipe.Scrr(:,j,k));
%         pipe2.Sctt(j,k) = mean(pipe.Sctt(:,j,k));
%         pipe2.Sczz(j,k) = mean(pipe.Sczz(:,j,k));
%         pipe2.Scrt(j,k) = mean(pipe.Scrt(:,j,k));
%         pipe2.Scrz(j,k) = mean(pipe.Scrz(:,j,k));
%         pipe2.Sctz(j,k) = mean(pipe.Sctz(:,j,k));
%               
        %Skewness tensor. Tensor of Rank 3. 
        
        pipe2.ururur(j,k) = mean(pipe.ururur(:,j,k));
        pipe2.ututut(j,k) = mean(pipe.ututut(:,j,k));
        pipe2.uzuzuz(j,k) = mean(pipe.uzuzuz(:,j,k));
        pipe2.ururut(j,k) = mean(pipe.ururut(:,j,k));
        pipe2.ururuz(j,k) = mean(pipe.ururuz(:,j,k));
        pipe2.urutut(j,k) = mean(pipe.urutut(:,j,k));
        pipe2.ututuz(j,k) = mean(pipe.ututuz(:,j,k));
        pipe2.uruzuz(j,k) = mean(pipe.uruzuz(:,j,k));
        pipe2.utuzuz(j,k) = mean(pipe.utuzuz(:,j,k));
        pipe2.urutuz(j,k) = mean(pipe.urutuz(:,j,k));

    end    
end

%Compute TKE budget terms
P_k  = 0.5*(pipe2.Prr  + pipe2.Ptt  + pipe2.Pzz);   %production
T_k  = 0.5*(pipe2.Trr  + pipe2.Ttt  + pipe2.Tzz);   %turbulent transport
PS_k = 0.5*(pipe2.PSrr + pipe2.PStt + pipe2.PSzz);  %pressure-strain
PT_k = 0.5*(pipe2.PTrr + pipe2.PTtt + pipe2.PTzz);  %pressure-transport
Pi_k =     (PT_k - PS_k);
VD_k = 0.5*(pipe2.VDrr + pipe2.VDtt + pipe2.VDzz);  %viscous diffusion
D_k  = 0.5*(pipe2.Drr  + pipe2.Dtt  + pipe2.Dzz);   %dissipation
C_k  = 0.5*(pipe2.Crr  + pipe2.Ctt  + pipe2.Czz);   %mean convection
K_k  = 0.5*(pipe2.urur + pipe2.utut + pipe2.uzuz)./(u_tau.^2);

%%Compute uTau
u_tau=sqrt(-pipe2.dUzdr(end,:)*nu);

%Write postprocessed results in a file
% make an array of radii
r_=zeros(nr,1);
for i=1:nr
    r_(i,1)=sqrt(x(i,1)^2+y(i,1)^2);
end
rMax=max(r_);
r_=rMax-r_;  %Changing the discrete radius from wall toward center

Re_tau=u_tau*rMax/nu;

for k=1:nz
yplus(:,k) = (Rmax-pipe.r1(1,:,k))*u_tau(k)/nu;
end

%Plot TKE
load 'Retau_180.dat'
rplus = Retau_180(:,2);
urms = Retau_180(:,5);
vrms = Retau_180(:,6);
wrms = Retau_180(:,7);
KK = 0.5*(urms.^2+vrms.^2+wrms.^2);
indices = [1 11 21 27 34 41 47 54 60 67 80 93];
figure()
for l = 1:length(indices)
    semilogy(K_k(:,indices(l))/4+zn(indices(l)),yplus(:,indices(l)),KK/4+zn(indices(l)),rplus,'k--')
    hold on
end
xlabel('z/D')
ylabel('(R-r)^+')
title('Turbulent Kinetic energy')

Ucl = max(pipe2.Uz);
Ucl_ref = 1.306579;
figure()
plot(zn,Ucl/Ucl_ref,zn,ones(length(zn),1),'r--','LineWidth',2)
title('Centerline velocity comparison')
hold on
shade(zn,0.99*ones(length(zn),1),'k',zn,1.01*ones(length(zn),1),'k','FillType',[1 2;2 1])

utau_ref = 6.8322301823104711e-02;
figure()
plot(zn, u_tau/utau_ref,zn,ones(length(zn),1),'r--','LineWidth',2)
title('Friction velocity comparison')
hold on
shade(zn,0.99*ones(length(zn),1),'k',zn,1.01*ones(length(zn),1),'k','FillType',[1 2;2 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% Shape factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% The shape factor in boundary layer flows is defined as:
%%% H12 = δ1/δ2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% δ1: is the displacement thickness  
%%% δ1 = int(0,inf)[1 - u(y)/uo]dy 
%%% uo is the velocity in the free-stream outside the boundary layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% δ2: is the momentum thickness 
%%% δ2 = int(0,inf)[u(y)/uo*{1 - u(y)/uo}]dy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Id1 = zeros(nz,1);
Id2 = zeros(nz,1);
delta1 = zeros(nz,1);
delta2 = zeros(nz,1);

for k = 1:nz
     Id1(k) = 2*trapz(pipe.r1(1,:,1)',pipe.r1(1,:,1)'.*(1-pipe2.Uz(:,k)/Ucl(k)));
     a1 = 1; b1 = -2*Rmax; c1 = Id1(k);
     delta1(k) = (-b1  - sqrt(b1^2 - 4*a1*c1))/(2*a1);
 
     Id2(k) = 2*trapz(pipe.r1(1,:,1)',pipe.r1(1,:,1)'.*(pipe2.Uz(:,k)/Ucl(k)).*(1-pipe2.Uz(:,k)/Ucl(k)));
     a2 = 1; b2 = -2*Rmax; c2 = Id2(k);
     delta2(k) = (-b2  - sqrt(b2^2 - 4*a2*c2))/(2*a2);
end

H12 = delta1./delta2;
H12_ref = 1.85;
figure()
plot(zn, H12/H12_ref,zn,ones(length(zn),1),'r--','LineWidth',2)
title('Shape factor comparison')
hold on
shade(zn,0.99*ones(length(zn),1),'k',zn,1.01*ones(length(zn),1),'k','FillType',[1 2;2 1])

% Mean velocity and rms plots
m = 70;
figure()
semilogx(yplus(:,m), pipe2.Uz(:,m)./u_tau(m))
title('Mean streamwise velocity')
ylabel('U_z^+')
xlabel('(R-r)^+')

figure()
semilogx(yplus(:,m), sqrt(pipe2.urur(:,m))/u_tau(m),...
     yplus(:,m), sqrt(pipe2.utut(:,m))/u_tau(m),...
     yplus(:,m), sqrt(pipe2.uzuz(:,m))/u_tau(m),...
     yplus(:,m), sqrt(pipe2.uruz(:,m))/u_tau(m))
legend('urur','utut','uzuz','uruz')
title('RMS values')
ylabel('urms^+')
xlabel('(R-r)^+')

figure()
plot(yplus(:,m), P_k(:,m)*nu/(u_tau(m)^4),...
    yplus(:,m), PT_k(:,m)*nu/(u_tau(m)^4),... 
    yplus(:,m), C_k(:,m)*nu/(u_tau(m)^4),...
        yplus(:,m), D_k(:,m)*nu/u_tau(m)^4,...
          yplus(:,m), VD_k(:,m)*nu/u_tau(m)^4)
    
legend('Production','Pressure-transport','Convection','Dissipation','Viscous Diffusion')
% legend('Production','Pressure-transport', 'Convection')
title('TKE Budget')
ylabel('k-budget')
xlabel('(R-r)^+')
%clearvars -except pipe3 yplus nu u_tau

% velocity profiles + derivatives
fOut=fopen('./stats_results/turbPipe_meanVelocity.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# Mean velocity and their derivatives \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# z\t uTau \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r\t z\t Ur\t Ut\t Uz\t dUrdr\t dUrdt\t dUrdz\t dUtdr\t dUtdt\t dUtdz\t dUzdr\t dUzdt\t dUzdz \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for j=1:nz
	for i = 1:nr
		fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t   \n',...
					r_(i,1),zn(j),pipe2.Ur(i,j),pipe2.Ut(i,j),pipe2.Uz(i,j),...
					pipe2.dUrdr(i,j),pipe2.dUrdt(i,j),pipe2.dUrdz(i,j),... 
					pipe2.dUtdr(i,j),pipe2.dUtdt(i,j),pipe2.dUtdz(i,j),... 
					pipe2.dUzdr(i,j),pipe2.dUzdt(i,j),pipe2.dUzdz(i,j))
	end						
end
fclose(fOut);

%Reynolds stress components
fOut=fopen('./stats_results/turbPipe_ReyStress.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# Reynolds stress components \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# z\t uTau \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r\t z\t urp\t urut\t uruz\t utp\t uzp\t utuz \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for j = 1:nz
     for i=1:nr
          fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t \n',...
                    r_(i,1),zn(j),sqrt(pipe2.urur(i,j)),pipe2.urut(i,j),pipe2.uruz(i,j),...
                    sqrt(pipe2.utut(i,j)),sqrt(pipe2.uzuz(i,j)),pipe2.utuz(i,j))
     end
end
fclose(fOut);

%TKE budget terms
fOut=fopen('./stats_results/turbPipe_kBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# TKE budget terms: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# z\t uTau \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t z+\t P_k+ \t T_k+ \t PS_k+ \t PT_k+ \t VD_k+ \t D_k+ \t C_k+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu./(u_tau.^4);
for j = 1:nz
     for i=1:nr
          fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t \n',...
                    r_(i,1)*u_tau(j)/nu,zn(j)*u_tau(j)/nu,P_k(i,j)*fac(j),T_k(i,j)*fac(j),PS_k(i,j)*fac(j),PT_k(i,j)*fac(j),VD_k(i,j)*fac(j),D_k(i,j)*fac(j),C_k(i,j)*fac(j))
     end
end
fclose(fOut);

%budget terms of <urur>
fOut=fopen('./stats_results/turbPipe_rrBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <urur>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# z\t uTau \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t z+\t P_rr+ \t T_rr+ \t PS_rr+ \t PT_rr+ \t VD_rr+ \t D_rr+ \t C_rr+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu./(u_tau.^4);
P_rr  = pipe2.Prr;   %production
T_rr  = pipe2.Trr;   %turbulent transport
PS_rr = pipe2.PSrr;  %pressure-strain
PT_rr = pipe2.PTrr;  %pressure-transport
VD_rr = pipe2.VDrr;  %viscous diffusion
D_rr  = pipe2.Drr;   %dissipation
C_rr  = pipe2.Crr;   %mean convection
for j=1:nz
     for i=1:nr
          fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t \n',...
                    r_(i,1)*u_tau(j)/nu,zn(j)*u_tau(j)/nu,P_rr(i,j)*fac(j),T_rr(i,j)*fac(j),PS_rr(i,j)*fac(j),PT_rr(i,j)*fac,VD_rr(i,j)*fac,D_rr(i,j)*fac(j),C_rr(i,j)*fac(j))
     end
end
fclose(fOut);

%budget terms of <utut>
fOut=fopen('./stats_results/turbPipe_ttBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <utut>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# z\t uTau \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t z+\t P_tt+ \t T_tt+ \t PS_tt+ \t PT_tt+ \t VD_tt+ \t D_tt+ \t C_tt+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu./(u_tau.^4);
P_tt  = pipe2.Ptt;   %production
T_tt  = pipe2.Ttt;   %turbulent transport
PS_tt = pipe2.PStt;  %pressure-strain
PT_tt = pipe2.PTtt;  %pressure-transport
VD_tt = pipe2.VDtt;  %viscous diffusion
D_tt  = pipe2.Dtt;   %dissipation
C_tt  = pipe2.Ctt;   %mean convection
for j=1:nz
     for i=1:nr
          fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t \n',...
                    r_(i,1)*u_tau(j)/nu,zn(j)*u_tau(j)/nu,P_tt(i,j)*fac(j),T_tt(i,j)*fac(j),PS_tt(i,j)*fac(j),PT_tt(i,j)*fac,VD_tt(i,j)*fac,D_tt(i,j)*fac(j),C_tt(i,j)*fac(j))
     end
end
fclose(fOut);

%budget terms of <uzuz>
fOut=fopen('./stats_results/turbPipe_zzBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <uzuz>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# z\t uTau \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t z+\t P_zz+ \t T_zz+ \t PS_zz+ \t PT_zz+ \t VD_zz+ \t D_zz+ \t C_zz+  \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu./(u_tau.^4);
P_zz  = pipe2.Pzz;   %production
T_zz  = pipe2.Tzz;   %turbulent transport
PS_zz = pipe2.PSzz;  %pressure-strain
PT_zz = pipe2.PTzz;  %pressure-transport
VD_zz = pipe2.VDzz;  %viscous diffusion
D_zz  = pipe2.Dzz;   %dissipation
C_zz  = pipe2.Czz;   %mean convection
for j=1:nz
     for i=1:nr
          fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t \n',...
                    r_(i,1)*u_tau(j)/nu,zn(j)*u_tau(j)/nu,P_zz(i,j)*fac(j),T_zz(i,j)*fac(j),PS_zz(i,j)*fac(j),PT_zz(i,j)*fac,VD_zz(i,j)*fac,D_zz(i,j)*fac(j),C_zz(i,j)*fac(j))
     end
end
fclose(fOut);

%budget terms of <uruz>
fOut=fopen('./stats_results/turbPipe_rzBudget.dat','w');
fprintf(fOut,'# Postprocessed results of turbulent pipe flow \n')
fprintf(fOut,'# budget terms of <uruz>: nondim by multiplying by nu/uTau^4 \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
for k = 1:nz
	fprintf(fOut,'%+4.8f\t %+4.8f\t \n', zn(k), u_tau(k))
end
fprintf(fOut,'# nu = %+4.8f \n',nu)
fprintf(fOut,'# ------------------------------------------------------------------\n')
fprintf(fOut,'# r+\t z+\t P_rz+ \t T_rz+ \t PS_rz+ \t PT_rz+ \t VD_rz+ \t D_rz+ \t C_rz+ \n')
fprintf(fOut,'# ------------------------------------------------------------------\n')
fac=nu./(u_tau.^4);
P_rz  = pipe2.Prz;   %production
T_rz  = pipe2.Trz;   %turbulent transport
PS_rz = pipe2.PSrz;  %pressure-strain
PT_rz = pipe2.PTrz;  %pressure-transport
VD_rz = pipe2.VDrz;  %viscous diffusion
D_rz  = pipe2.Drz;   %dissipation
C_rz  = pipe2.Crz;   %mean convection
for j=1:nz
     for i=1:nr
          fprintf(fOut,'%+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t %+4.8f\t \n',...
                    r_(i,1)*u_tau(j)/nu,zn(j)*u_tau(j)/nu,P_rz(i,j)*fac(j),T_rz(i,j)*fac(j),PS_rz(i,j)*fac(j),PT_rz(i,j)*fac,VD_rz(i,j)*fac,D_rz(i,j)*fac(j),C_rz(i,j)*fac(j))
     end
end
fclose(fOut);
