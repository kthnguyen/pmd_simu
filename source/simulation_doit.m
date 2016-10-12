%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation of TOF system using PMD  %
%%%% Author: Kien Nguyen                 %
%%%% Date: 18 Sep 16                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Synthesizing data %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constant
c = 299792458; %speed of light

%% Setup coordinates of components (in meter)

%PMD sensor
R_ssize = 45e-6;    %pixel dimension
R_sw = 160;         %num of width pixels
R_sh = 120;         %num of height pixels

R_x = 0;
R_y = -R_sw/2*R_ssize:R_ssize:R_sw/2*R_ssize-R_ssize;           %Center width at 0 m
R_z = (0.6-R_sh/2*R_ssize):R_ssize:(0.6+R_sh/2*R_ssize)-R_ssize;%Center height at 0.6 m

%Light source (single source with global illumination)
S_x = 0;
S_y = 1;
S_z = 0.6;

%Object (scene) which is observable by the PMD
%160x120 patches corresponding to 16x120 sensor pixels
O_psize = 1e-2; %patch dimesion

%rear panel
O_pw = 160;     %num of width patches
O_ph = 120;     %num of height patches
O_xr = 3;
O_yr = -O_pw/2*O_psize:O_psize:O_pw/2*O_psize-O_psize;
O_zr = 0:O_psize:O_psize*O_ph-O_psize;

%front panel
O_yflo = 41;
O_yfhi = 120;
O_zflo = 1;
O_zfhi = 60;

O_xf = 2.5;
O_yf = O_yr(O_yflo:O_yfhi);
O_zf = O_zr(O_zflo:O_zfhi);

%% Visualisation of the setup
%Rear panel
x = ones(R_sh,R_sw)*O_xr;
y = repmat(O_yr,R_sh,1);
z = repmat(O_zr,R_sw,1)';

figure(2)
surf(x,y,z,'LineStyle','none','FaceColor','k')
axis([-1 4 -1 1 0 2])
hold on

%Front panel
x = ones(60,80)*O_xf;
y = repmat(O_yf,60,1);
z = repmat(O_zf,80,1)';

surf(x,y,z,'LineStyle','none','FaceColor','y')
hold on

%PMD sensor, overscale for visualisation
Rsurf_ssize = 1e-3;
Rsurf_x = 0;
Rsurf_y = -R_sw/2*Rsurf_ssize:Rsurf_ssize:R_sw/2*Rsurf_ssize-Rsurf_ssize;           %Center width at 0 m
Rsurf_z = (0.6-R_sh/2*Rsurf_ssize):Rsurf_ssize:(0.6+R_sh/2*Rsurf_ssize)-Rsurf_ssize;%Center height at 0.6 m

x = ones(R_sh,R_sw)*Rsurf_x;
y = repmat(Rsurf_y,R_sh,1);
z = repmat(Rsurf_z,R_sw,1)';
surf(x,y,z,'LineStyle','none','FaceColor','r')

%Light source
plot3(S_x,S_y,S_z,'go','MarkerFaceColor','g','MarkerSize',10)


%% Calculate in-shade patches on the rear panel based on coordinates

%Rear panel plane function based on 3 patches
Or1 = [O_xr,O_yr(1),O_zr(1)];
Or2 = [O_xr,O_yr(end),O_zr(end)];
Or3 = [O_xr,O_yr(1),O_zr(end)];

O_normalr = cross(Or1-Or2,Or1-Or3);
syms x y z
O = [x,y,z];
O_planer = dot(O_normalr, O-Or1);%plane fn interms of x,y,z

%Line connecting the light source to 2 top corners patches of the front
Of1 = [O_xf, O_yf(1),O_zf(end)];
Of2 = [O_xf, O_yf(end),O_zf(end)];
Ss = [S_x, S_y, S_z];

syms t
O_linef1 = Ss + t*(Ss-Of1);%line fn in terms of t
O_linef2 = Ss + t*(Ss-Of2);

%Calculate intersection points at the rear panel
t1 = solve(subs(O_planer, O, O_linef1));%subs t into plane fn then solve
t2 = solve(subs(O_planer, O, O_linef2));

point1 = subs(O_linef1, t, t1);%x,y,z
point2 = subs(O_linef2, t, t2);%x,y,z

%indices of patches in shade so it can be called from O_yr and O_zr
sh_ylo = double(floor((point1(2)-O_yr(1))/O_psize)+1);
sh_yhi = double(floor((point2(2)-O_yr(1))/O_psize)+1);

sh_zlo = 1;
sh_zhi = double(floor(point1(3)/O_psize)+1);


%% Calculate time of flight for each patches observed by the PMD
%PMD is placed so that it is head straight to the scene. So the front panel
%occludes all rear panel pathces with same y and z coordinates

%Source to scene
%L dim: 160x120
L = sqrt(CalculateMutualSum((O_yr-S_y).^2,(O_zr-S_z).^2) + (O_xr-S_x)^2);%rear panel
L(O_yflo:O_yfhi,O_zflo:O_zfhi) = sqrt(CalculateMutualSum((O_yf-S_y).^2,(O_zf-S_z).^2) + (O_xf-S_x)^2);%front panel

%Scence to sensor
%D dim: 160x120
R_z = R_z(end:-1:1); %The image is inversed vertically through camera lens
R_y = R_y(end:-1:1);
D = sqrt(CalculateMutualSum((O_yr-R_y).^2,(O_zr-R_z).^2) + (O_xr-R_x)^2);%rear panel
D(O_yflo:O_yfhi,O_zflo:O_zfhi) = sqrt(CalculateMutualSum((O_yf-R_y(O_yflo:O_yfhi)).^2,(O_zf-R_z(O_zflo:O_zfhi)).^2) + (O_xf-R_x)^2);%front panel

%Total distance
K = L + D;

% Patches in shade gave 0 path
% DEAL WITH THE SITUATION WHEN PART OF THE SHADE IS COVERED BY THE FRONT
% PANEL IN THE SCENE OBSERVED BY THE PMD
% !(Need bugs checked for other positions of the light source)
if (sh_yhi < O_yflo)||(sh_ylo > O_yfhi) %not covered
    K(sh_ylo:sh_yhi,sh_zlo:sh_zhi) = 0;
elseif (sh_yhi > O_yflo)
    K(sh_ylo:O_yflo-1,sh_zlo:sh_zhi) = 0;
elseif (sh_ylo < O_yfhi)
    K(O_yfhi+1:sh_yhi,sh_zlo:sh_zhi) = 0;
end

% Time of flight 
Tau = K/c; %this makes up alpha(tau)
Tau = Tau(end:-1:1,:); %Inverse along width due to nature of camera system,height is ok

%% Construct and verify the expected transient image based on synthesized Tau
trans_img = zeros(160,120,600);

%max Tau
tau_max = max(Tau(:));
%min Tau due to some zero elements
tmp = Tau;
tmp(tmp==0) = Inf;
tau_min = min(tmp(:));

num_step = 600;
tau_step = (tau_max - tau_min)/(num_step-1);

trans_step = floor((Tau-tau_min)./tau_step)+1;
trans_step(trans_step<0) = 0;


% For each sensor pixel, its signal is 1 at corresponding time step in
% trans_step, then convolute with Gaussian function
w = gausswin(600,10);%narrow Guassian window
for m = 1:160
   for n = 1:120
       tmp = trans_step(m,n);
       if tmp > 0 % to avoid the shade patches
        trans_img(m,n,tmp) = 1;
       end
       trans_img(m,n,:) = conv(w,permute(trans_img(m,n,:),[3 1 2]),'same');
   end
end


% some rotations for displaying with image() fn
norm_img = permute(trans_img,[2 1 3]);
norm_img = norm_img(end:-1:1,:,:);
%norm_img(norm_img < 0.2) = 0.2; %ambient light


%% Show the transient image
while (1)
    for m = 1:600
        pause(3*tau_step);
        figure(10)
        image(norm_img(:,:,m)*200);
        colormap(gray(256));
    end
end

%% Set up the movie.
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
open(writerObj); 
 
for m = 1:600    
    pause(3*tau_step);
    figure(10)
    image(norm_img(:,:,m)*200);
    colormap(gray(256));
 
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
 
end
hold off
close(writerObj); % Saves the movie.
%%
figure (8)
stem(1:1:600,permute(norm_img(80,40,:),[3 1 2]))
title('Signal at pixel (80,40)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Implement Lin's algorithm %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%% With simple sinusoidal source %%%%%%%%%%%%%%%%%%%%
% Synthesizing H(w,phi)
Nt = 10; %Capture time: 10 signal periods
Amp = 1; %Source signal amplitude
fmin = 1; %MHz
fmax = 180; %MHz
fsource = (fmin:1:fmax)*1e6; %Source signal frequency range
fstep = 1;
f_len = (fmax - fmin)/fstep+1;
Bf = Nt*Amp^2./(2*fsource);
%Bf = 1;

tau_len = 350;
tau_step = 0.1; %ns

dataH = zeros(160,120,f_len,2);%pixel width, height; frequencies; 2 phases 0 and pi/2
for m = 1:160
    for n = 1:120
        temp_tau = Tau(m,n);
        dataH(m,n,:,1) =  cos(2*pi*fsource*temp_tau).*Bf;
        dataH(m,n,:,2) =  cos(2*pi*fsource*temp_tau + pi/2).*Bf; 
    end
end

%% Processing with Lin's algorithm
Hc = zeros(160,120,f_len);
Hc = dataH(:,:,:,1)+ 1i*dataH(:,:,:,2);

wHc = zeros(160,120,f_len);
trans_img_r = zeros(160,120,tau_len);

%window
beta = 6;
win = kaiser(f_len*2+1,beta);
win = win(end-f_len+1:end);

for m=1:f_len
    wHc(:,:,m) = Hc(:,:,m)*win(m)/Bf(m);
end

trans_img_r = ProcessFreqDomain(wHc,tau_step, tau_len);

%% Eliminate patches in shade
for m = 1:160
    for n = 1:120
        if(Tau(m,n)==0)
            trans_img_r(m,n,:) = zeros(1,1,tau_len);
        end
    end
end

%% Rotate the transient image for display with MATLAB function
norm_img_r = permute(trans_img_r,[2 1 3]);
norm_img_r = norm_img_r(end:-1:1,:,:);

%% Show transient image
while(1)
    for m = 100:tau_len
        pause(0.2*tau_step);
        figure(10)
        image(norm_img_r(:,:,m)*1400);
        colormap(gray(256));
    end
end

%% Test at some patches
%earliest
x = 45;%w
y = 20;%h
tr1 = permute(trans_img_r(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)

%% last
x = 160;%w
y = 2;%h
tr1 = permute(trans_img_r(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)

%% shade
x = 140;
y = 80;
tr1 = permute(norm_img_r(y,x,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)

%%
%%%%%%%%%%%%%%% With harmonics source %%%%%%%%%%%%%%%%%%%%
%% Load calibration matrix
cf_name = 'calib/PmCF6_stripe';
load(cf_name,'A1_str','phi1_str','A3_str','phi3_str','shutter');
%   A1, phi1: Amplitudes and phases of fundamental component
%             dimensions: frequency,0/90deg,w,h,shutter.
%   shutter£º shutter time.

s = length(shutter);
A1 = A1_str(:,:,1:160,:,s);
phi1 = phi1_str(:,:,1:160,:,s);
A3 = A3_str(:,:,1:160,:,s);
phi3 = phi3_str(:,:,1:160,:,s);

Ac1 = A1.*exp(1i*phi1);
Ac3 = A3.*exp(1i*phi3);

Bc1 = Ac1(:,1,:,:)+ Ac1(:,2,:,:);
Bc3 = Ac3(:,1,:,:)+ Ac3(:,2,:,:);

Bc1 = permute(Bc1,[1 3 4 2]);
Bc3 = permute(Bc3,[1 3 4 2]);
%% Frequency
fmin = 1; %MHz
fmax = 180; %MHz
fsource = (fmin:1:fmax)*1e6; %Source signal frequency range
fstep = 1;
f_len = (fmax - fmin)/fstep+1;
fsource = fsource(:);

tau_len = 800;
tau_step = 0.05; %ns

%% Synthesize H based on calibration matrix
H_0 = zeros(160,120,180);
H_90 = zeros(160,120,180);
Hc = zeros(160,120,180);
R = zeros(160,120,180);
%  m = 45;
%  n = 20;
for m = 1: 160
     for n = 1:120
        temp_t = Tau(m,n);
        Ac1_0 = Ac1(:,1,m,n);
        Ac3_0 = Ac3(:,1,m,n);
        Ac1_90 = Ac1(:,2,m,n);
        Ac3_90 = Ac3(:,2,m,n);
        %B1 = Ac1_0 + Ac1_90;
        %B1 = permute(B1(:,m,n),[2 3 1 4]);
       
        B1 = Bc1(:,m,n);
        B3 = Bc3(:,m,n);
        H_0(m,n,:) = Ac1_0.*exp(-1i*fsource*temp_t) + Ac3_0 .*exp(-1i*3*fsource*temp_t);         
        H_90(m,n,:) = Ac1_90.*exp(-1i*(fsource*temp_t+pi/2)) + Ac3_90.*exp(-1i*3*(fsource*temp_t+pi/2));
        Hc(m,n,:) = H_0(m,n,:)+1i*H_90(m,n,:);

        %Hc(m,n,:) = B1.*exp(-1i*fsource*2*pi*temp_t) + B3.*exp(-3*1i*fsource*2*pi*temp_t);
        %B1 = permute(Bc1(:,m,n),[2 3 1]);
        B1 = permute(Bc1(:,m,n),[2 3 1]);
        R(m,n,:) = Hc(m,n,:)./B1;
     end
 end

%% Rotate for windowing processing 
Bc1 = permute(Bc1,[2 3 1]);
Bc3 = permute(Bc3,[2 3 1]);

%%
beta = 15;
%win = kaiser(f_len*2+1,beta);
%win = win(end-f_len+1:end);
win = hamming(f_len);

for m=1:f_len
    wHc(:,:,m) = R(:,:,m)*win(m);
%     wBc1(:,:,m) = Bc1(:,:,m)*win(m);
%     wBc3(:,:,m) = Bc3(:,:,m)*win(m);
end

trans_img_r = ProcessFreqDomain(wHc,tau_step, tau_len);
%%
Bc31 = Bc3./(3*Bc1);
Bc31_t = ProcessFreqDomain(Bc31,tau_step, tau_len);
%% Eliminate patches in shade
for m = 1:160
    for n = 1:120
        if(Tau(m,n)==0)
            trans_img_r(m,n,:) = zeros(1,1,tau_len);
        end
    end
end
%% Rotate the transient image for display with MATLAB function
norm_img_r = permute(trans_img_r,[2 1 3]);
norm_img_r = norm_img_r(end:-1:1,:,:);

%% Show transient image without dilated removal
while(1)
    for m = 10:280
        pause(0.5*tau_step);
        figure(10)
        image(norm_img_r(:,:,m)*1200);
        colormap(gray(256));
    end
end

%% Dilated removal
trans_img_r_prc = zeros(size(trans_img_r));
%  m = 45;
%  n = 20;
for m = 1:160
     for n = 1:120
        pix = trans_img_r(m,n,:);
        pix_itp = interp(pix(:),3);
        Bc31_te = permute(Bc31_t(m,n,:),[3 2 1]);
        pix_itp_scl = conv(pix_itp(1:800),Bc31_te,'same');
        
        [pks,locs] = findpeaks(pix(:));
        locs = floor(locs/1.5);
        pix_to_rmv = zeros(800,1);
        pix_to_rmv(3*locs:3*locs*2) = pix_itp_scl(3*locs:3*locs*2);
        pix_prc = pix(:) - pix_to_rmv; 
%         pix_prc = pix(:) - pix_itp_scl; 
%         pix_prc(pix_prc>-0.08)= pix_prc(pix_prc>-0.08)+0.08;
%         pix_prc(pix_prc<-0.08)= 0;
        pix_prc(pix_prc<-0)= 0;
%         pix_prc_itp = interp(pix_prc(:),6);
        trans_img_r_prc(m,n,:) = permute(pix_prc,[3 2 1]);
     end
end
%% Eliminate patches in shade
for m = 1:160
    for n = 1:120
        if(Tau(m,n)==0)
            trans_img_r_prc(m,n,:) = zeros(1,1,tau_len);
        end
    end
end

%% Rotate the transient image for display with MATLAB function
norm_img_r_prc = permute(trans_img_r_prc,[2 1 3]);
norm_img_r_prc = norm_img_r_prc(end:-1:1,:,:);
%% Show transient image with dilated removal
while(1)
    for m = 10:200
        pause(1.5*tau_step);
        figure(10)
        image(norm_img_r_prc(:,:,m)*1200);
        colormap(gray(256));
    end
end


%% Compare non-processed and processed image
x = 160;%w
y = 2;%h
tr1 = permute(trans_img_r(x,y,:),[3 1 2]);
tr2 = permute(trans_img_r_prc(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
plot(t,tr1,'linewidth',2);
hold on
plot(t,tr2,'g','linewidth',2);
%set(p,'linewidth',2)

%% Test at some patches
%earliest
x = 45;%w
y = 20;%h
tr1 = permute(trans_img_r(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)
%% last
x = 160;%w
y = 2;%h
tr1 = permute(trans_img_r(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)
%% earliest
x = 45;%w
y = 20;%h
tr1 = permute(trans_img_r_prc(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)
%% last
x = 160;%w
y = 2;%h
tr1 = permute(trans_img_r_prc(x,y,:),[3 1 2]);
t = (0:tau_len-1)*tau_step;
figure; 
p = plot(t,tr1);
set(p,'linewidth',2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% THE POLARIMETER %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign Muller matrix
Mu = zeros(160,120,4,4);

Mx = [1 0 0 0;0 0.8 0 0;0 0 0.84 0; 0 0 0 0.6];
My = [1 0 0 0;0 0.8 0 0;0 0 0.87 0; 0 0 0 0.65];

M_temp = repmat(Mx,[1 1 160 120]);
Mu = permute(M_temp,[3 4 2 1]);

% sticker on front panel
for m = 60:100
    for n = 40:60
        Mu(m,n,:,:) = My;
    end
end

% sticker on rear panel
for m = 60:100
    for n = 100:120
        Mu(m,n,:,:) = My;
    end
end

%% Stokes matrices
S1 = [1;0;0;1];
T1 = [1;0;0;-1];
S2 = [1;0;0;1];
T2 = [1;0;0;1];

%% For S1,T1
% Synthesizing H(w,phi)
Nt = 10; %Capture time: 10 signal periods
Amp = 1; %Source signal amplitude
fmin = 1; %MHz
fmax = 180; %MHz
fsource = (fmin:1:fmax)*1e6; %Source signal frequency range
fstep = 1;
f_len = (fmax - fmin)/fstep+1;
Bf = Nt*Amp^2./(2*fsource);
%Bf = 1;

tau_len = 350;
tau_step = 0.1; %ns

dataH1 = zeros(160,120,f_len,2);%pixel width, height; frequencies; 2 phases 0 and pi/2
% m = 80;
% n = 55;
for m = 1:160
   for n = 1:120
        Mu_mn = reshape(Mu(m,n,:,:),[4 4]);
        I_factor = transpose(T1)*Mu_mn*S1;
        temp_tau = Tau(m,n);
        dataH1(m,n,:,1) =  cos(2*pi*fsource*temp_tau).*Bf*I_factor;
        dataH1(m,n,:,2) =  cos(2*pi*fsource*temp_tau + pi/2).*Bf*I_factor; 
   end
end

Hc1 = zeros(160,120,f_len);
Hc1 = dataH1(:,:,:,1)+ 1i*dataH1(:,:,:,2);

%Process
wHc1 = zeros(160,120,f_len);
trans_img_r1 = zeros(160,120,tau_len);

%window
beta = 6;
win = kaiser(f_len*2+1,beta);
win = win(end-f_len+1:end);

for m=1:f_len
    wHc1(:,:,m) = Hc1(:,:,m)*win(m)/Bf(m);
end

trans_img_r1 = ProcessFreqDomain(wHc1,tau_step, tau_len);

%% For S2,T2
% Synthesizing H(w,phi)
Nt = 10; %Capture time: 10 signal periods
Amp = 1; %Source signal amplitude
fmin = 1; %MHz
fmax = 180; %MHz
fsource = (fmin:1:fmax)*1e6; %Source signal frequency range
fstep = 1;
f_len = (fmax - fmin)/fstep+1;
Bf = Nt*Amp^2./(2*fsource);
%Bf = 1;

tau_len = 350;
tau_step = 0.1; %ns

dataH2 = zeros(160,120,f_len,2);%pixel width, height; frequencies; 2 phases 0 and pi/2
% m = 80;
% n = 55;
for m = 1:160
   for n = 1:120
        Mu_mn = reshape(Mu(m,n,:,:),[4 4]);
        I_factor = transpose(T2)*Mu_mn*S2;
        temp_tau = Tau(m,n);
        dataH2(m,n,:,1) =  cos(2*pi*fsource*temp_tau).*Bf*I_factor;
        dataH2(m,n,:,2) =  cos(2*pi*fsource*temp_tau + pi/2).*Bf*I_factor; 
   end
end

Hc2 = zeros(160,120,f_len);
Hc2 = dataH2(:,:,:,1)+ 1i*dataH2(:,:,:,2);

%Process
wHc2 = zeros(160,120,f_len);
trans_img_r2 = zeros(160,120,tau_len);

%window
beta = 6;
win = kaiser(f_len*2+1,beta);
win = win(end-f_len+1:end);

for m=1:f_len
    wHc2(:,:,m) = Hc2(:,:,m)*win(m)/Bf(m);
end

trans_img_r2 = ProcessFreqDomain(wHc2,tau_step, tau_len);

%% Combine Hc1 and Hc2 to get Hc_con
trans_img_r_con = zeros(size(trans_img_r2));
m = 40;
n = 25;
for m = 1:160
   for n = 1:120
    pk1 = max(trans_img_r1(m,n,:));
    pk2 = max(trans_img_r2(m,n,:));
    pk_con = abs((pk1-pk2)/(pk1+pk2));
    trans_img_r_con(m,n,:) = trans_img_r1(m,n,:)/pk1*pk_con;
   end
end

%% Eliminate patches in shade
for m = 1:160
    for n = 1:120
        if(Tau(m,n)==0)
            trans_img_r_con(m,n,:) = zeros(1,1,tau_len);
            trans_img_r1(m,n,:) = zeros(1,1,tau_len);
            trans_img_r2(m,n,:) = zeros(1,1,tau_len);
        end
    end
end

%% Rotate the transient image for display with MATLAB function
norm_img_r_con = permute(trans_img_r_con,[2 1 3]);
norm_img_r_con = norm_img_r_con(end:-1:1,:,:);

norm_img_r1 = permute(trans_img_r1,[2 1 3]);
norm_img_r1 = norm_img_r1(end:-1:1,:,:);

norm_img_r2 = permute(trans_img_r2,[2 1 3]);
norm_img_r2 = norm_img_r2(end:-1:1,:,:);

%% Show transient image
while(1)
    
    for m = 100:300
        pause(0.5*tau_step);
        fig2 = figure(9);
        image(norm_img_r2(:,:,m)*500);
        colormap(gray(256));
        set(fig2, 'Position', [100, 400, 160*3, 120*3])
    end
    
    
      for m = 100:300
        pause(0.5*tau_step);
        fig1 = figure(8);
        image(norm_img_r1(:,:,m)*1200);
        colormap(gray(256));
        set(fig1, 'Position', [700, 400, 160*3, 120*3])
      end
      
    for m = 100:300
        pause(0.5*tau_step);
        fig_con = figure(10);
        image(norm_img_r_con(:,:,m)*300);
        colormap(gray(256));
        set(fig_con, 'Position', [1300, 400, 160*3, 120*3])
    end
      
      
end

%% Show constrast image only
  for m = 100:250
        pause(0.5*tau_step);
        fig_con = figure(10);
        image(norm_img_r_con(:,:,m)*300);
        colormap(gray(256));
        set(fig_con, 'Position', [1000, 400, 160*3, 120*3])
    end
      
%% Show S1,T1
while(1)
    for m = 100:250
        pause(0.3*tau_step);
        fig1 = figure(8);
        image(norm_img_r1(:,:,m)*1000);
        colormap(gray(256));
        set(fig1, 'Position', [700, 400, 160*3, 120*3])
    end
end


%% Show S2,T2
while(1)
    for m = 100:250
        pause(0.5*tau_step);
        fig2 = figure(9);
        image(norm_img_r2(:,:,m)*500);
        colormap(gray(256));
        set(fig2, 'Position', [700, 400, 160*3, 120*3])
    end
end


%%
figure
% plot(pix(:))
plot(reshape(trans_img_r_con(40,25,:),[1 350]))
hold on
plot(reshape(trans_img_r_con(160,2,:),[1 350]),'r--')
%%
figure
% plot(pix(:))
plot(reshape(trans_img_r1(40,25,:),[1 350]))
hold on
plot(reshape(trans_img_r2(40,25,:),[1 350]),'r--')
plot(reshape(trans_img_r_con(40,25,:),[1 350]),'g--')
%% For reference what constrast image should looked like
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesizing H(w,phi)
Nt = 10; %Capture time: 10 signal periods
Amp = 1; %Source signal amplitude
fmin = 1; %MHz
fmax = 180; %MHz
fsource = (fmin:1:fmax)*1e6; %Source signal frequency range
fstep = 1;
f_len = (fmax - fmin)/fstep+1;
Bf = Nt*Amp^2./(2*fsource);
%Bf = 1;

tau_len = 350;
tau_step = 0.1; %ns

dataH_con = zeros(160,120,f_len,2);%pixel width, height; frequencies; 2 phases 0 and pi/2
% m = 80;
% n = 55;
for m = 1:160
   for n = 1:120
        Mu_mn = reshape(Mu(m,n,:,:),[4 4]);
        I1 = transpose(T1)*Mu_mn*S1;
        I2 = transpose(T2)*Mu_mn*S2;
        I_factor_con = abs((I1-I2)/(I1+I2));
        temp_tau = Tau(m,n);
        dataH_con(m,n,:,1) =  cos(2*pi*fsource*temp_tau).*Bf*I_factor_con;
        dataH_con(m,n,:,2) =  cos(2*pi*fsource*temp_tau + pi/2).*Bf*I_factor_con; 
   end
end

Hc_con = zeros(160,120,f_len);
Hc_con = dataH_con(:,:,:,1)+ 1i*dataH_con(:,:,:,2);

%% Constrast image
dataH_con = zeros(160,120,f_len,2);

Hc_con = zeros(160,120,f_len);
 m = 45;
 n = 20;
for m = 1: 160
    for n = 1:120
    H1_0 = dataH1(m,n,:,1);
    H1_90 = dataH1(m,n,:,2);
    H2_0 = dataH2(m,n,:,1);
    H2_90 = dataH2(m,n,:,2);
    dataH_con(m,n,:,1) = abs((H1_0-H2_0)./(H1_0+H2_0));
    dataH_con(m,n,:,2) = abs((H1_90-H2_90)./(H1_90+H2_90));
    end
end

Hc_con = dataH_con(:,:,:,1)+ 1i*dataH_con(:,:,:,2);


%% Processing with Lin's algorithm
% Hc = zeros(160,120,f_len);
% Hc = dataH(:,:,:,1)+ 1i*dataH(:,:,:,2);

wHc = zeros(160,120,f_len);
trans_img_r = zeros(160,120,tau_len);

%window
beta = 6;
win = kaiser(f_len*2+1,beta);
win = win(end-f_len+1:end);

for m=1:f_len
    wHc(:,:,m) = Hc_con(:,:,m)*win(m)/Bf(m);
end

trans_img_r = ProcessFreqDomain(wHc,tau_step, tau_len);

%% Eliminate patches in shade
for m = 1:160
    for n = 1:120
        if(Tau(m,n)==0)
            trans_img_r(m,n,:) = zeros(1,1,tau_len);
        end
    end
end

%% Rotate the transient image for display with MATLAB function
norm_img_r = permute(trans_img_r,[2 1 3]);
norm_img_r = norm_img_r(end:-1:1,:,:);

%% Show transient image
while(1)
    for m = 100:tau_len
        pause(0.5*tau_step);
        figure(10)
        image(norm_img_r(:,:,m)*1400);
        colormap(gray(256));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%
%pix = trans_img_r(40,25,:);
figure
% plot(pix(:))
plot(reshape(trans_img_r(40,25,:),[1 350]))
hold on
plot(reshape(trans_img_r(160,2,:),[1 350]),'r--')




%% %%%%%%%%%%%%%%%%%%% FOR TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
m = 160;
n = 2;
figure
plot(reshape(dataH1(m,n,:,1),[1 180]),'r')
hold on
plot(reshape(dataH2(m,n,:,1),[1 180]),'b--')

%%
figure
plot(reshape(Hc_con(40,25,:),[1 180]))
hold on
plot(reshape(Hc_con(160,2,:),[1 180]),'r--')
%%
figure
%plot(reshape(Hc_con(40,25,:),[1 180]))
hold on
plot(reshape(Hc1(40,25,:),[1 180]),'r')
plot(reshape(Hc2(40,25,:),[1 180]),'g--')


%%
trans_img_r_prc = zeros(size(trans_img_r));
 m = 45;
 n = 20;
% for m = 1:160
%      for n = 1:120
        pix = trans_img_r(m,n,:);
        pix_itp = interp(pix(:),3);
        Bc31_te = permute(Bc31_t(m,n,:),[3 2 1]);
        pix_itp_scl = conv(pix_itp(1:800),Bc31_te,'same');
        [pks,locs] = findpeaks(pix(:));
        locs = locs/1.5;
        pix_to_rmv = zeros(800,1);
        pix_to_rmv(3*locs:3*locs*2) = pix_itp_scl(3*locs:3*locs*2);
        pix_prc = pix(:) - pix_to_rmv; 
%         pix_prc(pix_prc>-0.08)= pix_prc(pix_prc>-0.08)+0.08;
%         pix_prc(pix_prc<-0.08)= 0;
        %pix_prc(pix_prc<-0)= 0;
%         pix_prc_itp = interp(pix_prc(:),6);
        trans_img_r_prc(m,n,:) = permute(pix_prc,[3 2 1]);
%      end
% end


%%
figure
plot(pix(:))
hold on
plot(pix_itp(1:800),'g')
plot( pix_itp_scl,'r')
plot(pix_prc,'k')


%%
x = trans_img_r(45,20,:);
y = interp(x(:),3);
Bc31_te = permute(Bc31_t(45,20,:),[3 2 1]);
z = conv(y(1:800),Bc31_te,'same');
%k = deconv(z,Bc1_te);

figure
plot(x(:))
hold on
plot(z,'r')
plot(y,'g')
plot(Bc31_te,'k')
%plot(alpha)
%plot(y(1:800))

figure
plot(x(:)-z)

%%
figure
plot(Bc1_te)

figure
plot(Bc3_te)

%%
figure
stem(1:1:600,permute(trans_img(80,40,:),[3 1 2]))

%%
figure
stem(1:1:600,permute(trans_img_r(80,40,:),[3 1 2]))
%%
figure
stem(permute(ifft(Hc(80,40,:)),[3 1 2]))

%%
w = gausswin(600,10);
out = conv(w,permute(trans_img(1,1,:),[3 1 2]),'same');
figure
stem(out)

%%
m = 80;
n = 40;
 temp_tau = Tau(m,n);
        dataH(m,n,:,1) =  cos(2*pi*fsource*temp_tau)*Nt*Amp./(2*fsource);


%% Reconstructing shade at Tau(140,40)
% Synthesizing H(w,phi)
Nt = 10; %Capture time: 10 signal periods
Amp = 1; %Source signal amplitude
fmin = 1; %MHz
fmax = 180; %MHz
fsource = (fmin:1:fmax)*1e6; %Source signal frequency range
%Bf = Nt*Amp^2./(2*fsource);
Bf=1;

dataH = zeros(160,120,180,2);%pixel width, height; frequencies; 2 phases 0 and pi/2
m = 140;
n = 40;
temp_tau = Tau(m,n)
dataH(m,n,:,1) =  cos(2*pi*fsource*temp_tau);
dataH(m,n,:,2) =  cos(2*pi*fsource*temp_tau + pi/2);


% Processing with Lin's algorithm
Hc = zeros(160,120,180);
Hc(m,n,:) = dataH(m,n,:,1)+ 1i*dataH(m,n,:,2);

trans_img_r = zeros(160,120,600);
temp = permute(Hc(m,n,:),[1 3 2]);
% temp = [temp];
% Bf = [Bf];
% win = hamming(180);

beta = 6;
 win = kaiser(180*2+1,beta);
 win = win(end-180+1:end);
%win = kaiser(180);
%f_signal = temp./Bf;
pix = temp.*win';

pix = pix(:);
A = [abs(pix(1));pix];
a = ifft(A,20000);
a = 2*abs(real(a(1:600)));
trans_img_r(m,n,:) = permute(a,[3 2 1]);

figure
plot(a)
%%
temp = dataH(m,n,:,2);
temp = temp(:);
figure
plot(abs(temp))
%% 
figure(2)
subplot 211
stem(1:1:600,permute(trans_img(80,40,:),[3 1 2]))

subplot 212
stem(permute(trans_img_r(80,40,:),[3 1 2]))


