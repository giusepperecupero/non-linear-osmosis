%% MATLAB Codes for the NON-LINEAR OSMOSIS - SHADOW/LIGHT REMOVAL
%  Copyright (c) 2022, Giuseppe Antonio Recupero
%  All rights reserved.
%
%  Author:
%  Giuseppe Antonio Recupero (email: giuseppe.recupero3@unibo.it)
%
%  Address:
%  Department of Mathematics
%  University of Bologna
%  Piazza di Porta San Donato 5
%  40126, Bologna, Italy
%
%  Date:
%  June, 2022

%% Related folders in \dataset_removal :
% input : shadowed images
% mask : shadow edges
% original : ground-truth images (used only if GT images are known)
% output : results

% \output folder includes subfolders for each input, with:
% input image, edge mask with adapted thickness, 
% output image named after all the parameters

%%
clear all, close all

% image input
name = 'desk6_8.png';
name_edge = 'desk6.png';

% parameters of evolutive scheme
tau = 1e3; 
% parameters of the model
eps = 1e-7; offset = 0.01; p = 1;
% edge thickness correction
erode_edge = 0;
dilate_edge = 6;
edge = - erode_edge + dilate_edge;

% model: linear or non-linear
flag_linear = true;
% evolution scheme: explicit or semi-implicit
semiimplicit = true; 

% stopping criterion
stop_flag = 3; % 0 for relative change, 1 for SSIM, 2 for MSE, 3 for max_IT
rel_ch_min = 1e-8;
SSIM_max = 0.955;
MSE_min = 7e-4;
time = [0,1e4]; % maxit = time(2)/tau
show_plots = false;

if stop_flag == 1 || stop_flag == 2 % if groundtruth is known
    name_orig = '4.png';
    u_orig = double(imread(strcat('GT/',name_orig)))./255;
end

cd = cd;
folder = strcat(cd(1:end-5),'dataset_removal');
mkdir(strcat(folder,'/output/',name(1:end-4)));
%% READ DATA + NORMALISATION
umat = double(imread(strcat(folder,'/input/',name)));
umat = umat./255;
%umat = rgb2gray(umat);
[mx, my, c] = size(umat);

imwrite(umat,strcat(folder,'/output/',name(1:end-4),'/',name(1:end-4),'_input','.png'),'png')

umat = umat+offset; % to avoid singularities

umask_edge = double(imread(strcat(folder,'/mask/',name_edge)))./255;
umask_edge = umask_edge(:,:,1);
umask_edge=double(~imerode(double(~umask_edge),ones(erode_edge+1)));
umask_edge=double(~imdilate(double(~umask_edge),ones(dilate_edge+1)));
imwrite(umask_edge,strcat(folder,'/output/',name(1:end-4),'/',name(1:end-4),'_size_',string(edge),'.png'),'png')

% Show input, mask and evolution of the osmosis process
figure(100)
subplot(1,3,2)
imshow(umat-offset)
subplot(1,3,1)
imshow(imfuse(umat-offset,umask_edge))

%% INITIALIZATION

% UNCOMMENT FOR OSMOSIS
%av = sum(sum(sum(umat)))/(mx*my*c);
%u_new = av*ones(mx,my,c);
%umask_edge = ones(size(umask_edge)); 

% UNCOMMENT FOR SHADOW REMOVAL
u_new = umat;

showfig = 1; % #iteration between successive imshow
showfig = showfig*c;
tic

%% Precompute the drift term d = \nabla v / v
hx = 1; hy = 1;
for k = 1:c
    d1ij = zeros(mx+1,my,c);
    d2ij = zeros(mx,my+1,c);
    d1ij(2:mx,:,:) = diff(umat,1,1)./(umat(2:mx,:,:)+umat(1:mx-1,:,:))*2/hx;
    d2ij(:,2:my,:) = diff(umat,1,2)./(umat(:,2:my,:)+umat(:,1:my-1,:))*2/hy;

    d1ij(padarray(umask_edge,[1,0],'replicate','pre')==0) = 0;
    d2ij(padarray(umask_edge,[0,1],'replicate','pre')==0) = 0;
end

%% ITERATIONS
it = 0;
rel_ch_list = [];
SSIM_list = [];
MSE_list =[];

for tautimesk = time(1):tau:time(end) % until maxit is reached
    u_old = u_new; 
    for k=1:c % for each color channel
        % compute matrix A
        A = divmatrix(mx,my,u_old(:,:,k),umat(:,:,k),d1ij(:,:,k),d2ij(:,:,k),eps,p,umask_edge,flag_linear);
        
        if semiimplicit % u_new = (I-tau*A)^-1 * u_old
            Q = speye(mx*my)-tau*A; 
            u_out = Q\reshape(u_old(:,:,k),mx*my,1);
            % u_out = bicgstab(Q,u_lin);
        else % u_new = (I+tau*A) * u_old
            Q = speye(mx*my)+tau*A;
            tau_sup = 1/max(diag(abs(A))) % print stability constraint
            u_out = Q*reshape(u_old(:,:,k),mx*my,1);
        end
        u_new(:,:,k) = reshape(u_out,mx,my); % update u_new on channel k
        
        it = it+1;
        if floor(it/showfig)== it/showfig % show updated u_new
            figure(100), subplot(1,3,3)
            imshow(u_new-offset,[0 1])
            title(['Time step t=' num2str(tautimesk)])
            %pause(0.1)
        end
    end
    
    % Stopping criteria
    switch stop_flag
        
        case 0 % relative change
            rel_ch = sum(sum(sum(abs(u_new-u_old))))/(mx*my*c);
            if show_plots
                rel_ch_list = [rel_ch_list,rel_ch];
                figure(10), plot(rel_ch_list), title('Relative change')
            end
            if rel_ch < rel_ch_min
                fprintf('\n Number Iterations: %0.4f\n', tautimesk/tau);
                toc
            return
            end
            
        case 1 % SSIM
            SSIM = ssim(u_new-offset,u_orig);
            if show_plots
                SSIM_list = [SSIM_list, SSIM];
                figure(22), plot(SSIM_list), title('SSIM error')
            end
            if SSIM > SSIM_max
                fprintf('\n Number Iterations: %0.2f\n', tautimesk/tau);
                toc
                return
            end
            
        case 2 % MSE
            MSE = immse(u_new-offset,u_orig);
            if show_plots
                MSE_list = [MSE_list,MSE]; 
                figure(23), plot(MSE_list), title('MSE error')
            end
            if MSE < MSE_min
                fprintf('\n Number Iterations: %0.2f\n', tautimesk/tau);
                toc
                return
            end
    end
end

U = u_new-offset;
imshow(U,[0 1])

if flag_linear
    model = 'lin';
else
    model = 'g';
end

if semiimplicit
    method = 'SI';
else
    method = 'I';
end

imwrite(U,strcat(folder,'/output/',name(1:end-4),'/',name(1:end-4),'_edge',string(edge),'_p',string(p),'_',model,method,'_tau',string(tau),'_T',string(tautimesk),'_eps',string(eps),'.png'),'png')