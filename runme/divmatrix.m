function A = divmatrix(mx,my,u,v,d1ij,d2ij,epsilon,p,umask_edge,flag_linear)

g = @(z,epsilon,p) 1./(z.^p+epsilon); % diffusivity function
hx = 1; hy = 1; c = 1; % space discretization step

%% discretization of g on each pixel

g_out = ones(mx,my); % for linear model

if ~flag_linear % for non-linear model
    
    g_vert = zeros(mx,my);
    g_vert(2:mx-1,:) = (u(1:mx-2,:)-u(3:mx,:))/(2*hy) ;
    secterm_vert = zeros(mx,my);
    secterm_vert(2:mx-1,:) = (v(1:mx-2,:)-v(3:mx,:)).*u(2:mx-1,:)./(2*hy*v(2:mx-1,:)) ;
    g_vert(umask_edge==1) = g_vert(umask_edge==1)-secterm_vert(umask_edge==1);

    g_hor = zeros(mx,my);
    g_hor(:,2:my-1) = (u(:,1:my-2)-u(:,3:my))/(2*hx) ;
    secterm_hor = zeros(mx,my);
    secterm_hor(:,2:my-1) = (v(:,1:my-2)-v(:,3:my)).*u(:,2:my-1)./(2*hx*v(:,2:my-1)) ;
    g_hor(umask_edge==1) = g_hor(umask_edge==1)-secterm_hor(umask_edge==1);

    g_discr = sqrt(g_vert.^2 + g_hor.^2); % || \nabla u - d*u ||
    g_out = g(g_discr,epsilon,p);
end

%% computation of A

% gradient upper (u_{i+1,j} - u_{ij})  and lower (u_{ij}-u_{i-1,j})
grad1xup  = spdiags([-ones(mx,1) ones(mx,1)],[0,1],mx,mx);
grad1xlow = spdiags([-ones(mx,1) ones(mx,1)],[-1,0],mx,mx);
% gradient upper (u_{i,j+1} - u_{ij}) and lower (u_{ij}-u_{i,j-1})
grad1yup  = spdiags([-ones(my,1) ones(my,1)],[0,1],my,my);
grad1ylow = spdiags([-ones(my,1) ones(my,1)],[-1,0],my,my);

Grad1xup  = kron(speye(my),grad1xup);
Grad1xlow = kron(speye(my),grad1xlow);
Grad1yup  = kron(grad1yup,speye(mx));
Grad1ylow = kron(grad1ylow,speye(mx));

d1ij(padarray(umask_edge,[1,0],'replicate','pre')==0) = 0;
d2ij(padarray(umask_edge,[0,1],'replicate','pre')==0) = 0;

g1ij = zeros(mx+1,my,c);
g2ij = zeros(mx,my+1,c);

g1ij(2:mx,:,:) = (g_out(2:mx,:,:)+g_out(1:mx-1,:,:))/2;
g2ij(:,2:my,:) = (g_out(:,2:my,:)+g_out(:,1:my-1,:))/2;

g_dx = spdiags(reshape(g1ij(2:mx+1,:,1),mx*my,1),0,mx*my,mx*my);
g_sx = spdiags(reshape(g1ij(1:mx,:,1),mx*my,1),0,mx*my,mx*my);
g_low = spdiags(reshape(g2ij(:,2:my+1,1),mx*my,1),0,mx*my,mx*my);
g_up = spdiags(reshape(g2ij(:,1:my,1),mx*my,1),0,mx*my,mx*my);

% average upper (u_{i+1,j} + u_{ij})/2  and lower (u_{ij}+u_{i-1,j})/2
m1xup  = spdiags(ones(mx,2)/2,[0,1],mx,mx);
m1xlow = spdiags(ones(mx,2)/2,[-1,0],mx,mx);
% average upper (u_{i,j+1} + u_{ij})/2  and lower (u_{ij}+u_{i,j-1})/2
m1yup  = spdiags(ones(my,2)/2,[0,1],my,my);
m1ylow = spdiags(ones(my,2)/2,[-1,0],my,my);

M1xup  = kron(speye(my),m1xup);
M1xlow = kron(speye(my),m1xlow);
M1yup  = kron(m1yup,speye(mx));
M1ylow = kron(m1ylow,speye(mx));

A  = 1/hx*( g_dx * ( Grad1xup - spdiags(reshape(d1ij(2:mx+1,:,1),mx*my,1),0,mx*my,mx*my)*M1xup ) -...
    g_sx * ( Grad1xlow - spdiags(reshape(d1ij(1:mx,:,1),mx*my,1),0,mx*my,mx*my)*M1xlow))+...
    +1/hy*( g_low * ( Grad1yup - spdiags(reshape(d2ij(:,2:my+1,1),mx*my,1),0,mx*my,mx*my)*M1yup )-...
    g_up * ( Grad1ylow - spdiags(reshape(d2ij(:,1:my,1),mx*my,1),0,mx*my,mx*my)*M1ylow));
end


