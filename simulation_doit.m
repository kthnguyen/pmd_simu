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
Tau = K/c;
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
trans_img(trans_img < 0.2) = 0.2; %ambient light

% some rotations for displaying with image() fn
trans_img = permute(trans_img,[2 1 3]);
trans_img = trans_img(end:-1:1,:,:);


%% Show the transient image
for m = 1:600
    pause(0.2*tau_step);
    figure(10)
    image(trans_img(:,:,m)*200);
    colormap(gray(256));
end
    
figure (8)
stem(1:1:600,permute(trans_img(80,40,:),[3 1 2]))
title('Signal at pixel (80,40)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Implement Lin's algorithm %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%% %%%%%%%%%%%%%%%%%%% FOR TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure
stem(1:1:600,permute(trans_img(80,40,:),[3 1 2]))

%%
w = gausswin(600,10);
out = conv(w,permute(trans_img(1,1,:),[3 1 2]),'same');
figure
stem(out)

%%
figure
stem(w)






