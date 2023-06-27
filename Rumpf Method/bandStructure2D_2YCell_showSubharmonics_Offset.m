
% bandStructure2 calculates the bandstructure for an infinite photonic
% crystal of unit cells defined by unitCell. 
% unitCell(1,:,:) sets the cell εᵣ, and unitCell(2,:,:) the μᵣ. 
% The number of harmonics used in the calculation are P and Q for the x and
% y dimensions respectively, and the number of points calculated along the 
% edges of the brillouin zone is set by N.

% This script in particular calculates the band structure of a
% non-primitive unit cell. I.e. a unit cell which contains more than 1
% primitive unit cells.
% In this case, the unit cell consists of two square 'air-hole' unit cells
% stacked on top of each other in the y axis.
% Further, this script offsets the points selected in reciprocal space, to
% show that the non-physical bands are actually the same bandstructure,
% just offset from that of the primitive cell.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup

startTime = tic();
pool = gcp();

% Define lattice vectors and reciprocal lattice vectors
a1_ = [1; 0; 0];
a2_ = [0; 1; 0];

Lx = norm(a1_);
Ly = norm(a2_);

P = 17;
Q = 9;
NH = P*Q; % total number of harmonics

N = 200;        % The number of bloch vectors to calculate between points in reciprocal space.

er = 8.9;       % The permittivity of the dielectric material
width = 0.835;  % The fraction of the unit cell covered by the square air hole

resX = 1000;    % The resolution of the unit cell, in X and Y.
resY = 1000;

closeFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add 'Utility Scripts' to path, if needed %%
if exist('Utility Scripts','dir') == 7
    % directory found, do nothing.
else
    fprintf('Finding ''Utility Scripts'' Folder\n');
    this_folder = cd;
    % loop to find parent directory containing 'Utility Scripts' folder
    searching = true;
    while searching
        current_folder = pwd;
        fprintf('path: %s\n',pwd);
        if isfolder('Utility Scripts')
            % folder found!
            cd 'Utility Scripts'
            addpath(genpath(pwd));
            cd(this_folder)
            searching = false;
        else
            cd ..
            if strcmp(current_folder, pwd)
                error('Error: Unable to locate ''Utility Scripts'' directory in any parent directory!');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define unit cell %%
eGrid = er*ones(resY,resX);
eGrid(floor(resY*(0.5-(width/(Ly*2)))):floor(resY*(0.5+(width/(Ly*2)))),floor(resX*(0.5-(width/(Lx*2)))):floor(resX*(0.5+(width/(Lx*2))))) = 1;

uGrid = ones(resY,resX);

% CHANGE UNIT CELL AND VALUES TO HAVE DOUBLE UNIT CELL!
eGrid = [eGrid(:,:);eGrid(:,:)];
uGrid = [uGrid(:,:);uGrid(:,:)];
a2_ = 2*a2_;

Lx = norm(a1_);     %recalculate lattice constants.
Ly = norm(a2_);
% END CHANGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute convolution matrices %%
fprintf('Computing Convolution Matrices');
startConv = tic();

Er = convMat(eGrid,P,Q);
Ur = convMat(uGrid,P,Q);

fprintf([': ', num2str(toc(startConv)), 's\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct list of points in the 1st Brillouin Zone to path between %%

% Compute reciprocal lattice vectors
ac = norm(cross(a1_,a2_));  % area of unit cell.

ahash_ = (2*pi/(ac)) * cross(a2_,[0;0;1]);  %create vector orth. to lattice vector
bhash_ = (2*pi/(ac)) * cross([0;0;1],a1_);

ahash_ = ahash_(1:2);   %reduce to 2D
bhash_ = bhash_(1:2);

% Construct list of BZ points
points = [{[0;0],"\Gamma"}; {ahash_/2,"X"}; {(ahash_+2*bhash_)/2,"M"}];

% Offset points by bhash_, to demonstrate the nonphysical bands are just
% the regular bands, from a different starting location.
points(:,1) = cellfun(@(x) x-bhash_,points(:,1),'un',0);
points(:,2) = cellfun(@(x) x+"'",points(:,2),'un',0);

[numPoints,~] = size(points);       % How many points we want to go between

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute list of betas %%

b = [];
dists = [0];
interval = (1/(2*N):1/N:1);
for i = 1:length(points)
    p1 = points(i,:);
    if i == length(points)
        p2 = points(1,:);
    else
        p2 = points(i+1,:);
    end
    
    delta = (p2{1}-p1{1})*(1/(2*N):1/N:1);
    p1top2_ = p1{1} + delta;
    
    dist = dists(end) + sqrt(diag(delta.' * delta)).';
    b = [b p1top2_];
    dists = [dists dist];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate PWE for each beta %%
fprintf('Calculating Betas');
startBand = tic();

central_index = (P+1)/2;
physical_indices = central_index + 2*(-floor(P/4):floor(P/4));

NSH_physical = (1+floor(P/4)*2)*Q;    % Number of subharmonics - periodic on double cell, but not primitive cell
NSH_nonphysical = NH-NSH_physical;

SHInds_phys = reshape(physical_indices'+(0:1:Q-1)*P,[1,NSH_physical]);
SHInds_nonphys = setdiff(1:P*Q,SHInds_phys);


AEs = zeros(numPoints*N,NH,NH);
AHs = zeros(numPoints*N,NH,NH);
BEs = zeros(numPoints*N,NH,NH);
BHs = zeros(numPoints*N,NH,NH);

VE_dat_p = zeros(numPoints*N,NSH_physical,NSH_physical);    % Arrays to store the eigenvectors
VH_dat_p = zeros(numPoints*N,NSH_physical,NSH_physical);
DE_dat_p = zeros(numPoints*N,NSH_physical);
DH_dat_p = zeros(numPoints*N,NSH_physical);

VE_dat_np = zeros(numPoints*N,NSH_nonphysical,NSH_nonphysical);
VH_dat_np = zeros(numPoints*N,NSH_nonphysical,NSH_nonphysical);
DE_dat_np = zeros(numPoints*N,NSH_nonphysical);
DH_dat_np = zeros(numPoints*N,NSH_nonphysical);


beta_xs = b(1,:);
beta_ys = b(2,:);

% Construct array of G-G' vectors
p = -((P-1)/2):+((P-1)/2); %indices along x
q = -((Q-1)/2):+((Q-1)/2); %indices along y
clear KXs KYs KX KY KXYs
KXYs(1,:,:) = zeros(P,Q);
KXYs(2,:,:) = zeros(P,Q);


for i = 1:P
    for j = 1:Q
        KXYs(:,i,j) = bhash_*p(i) + ahash_*q(j);
    end
end


KXs = squeeze(KXYs(1,:,:));
KYs = squeeze(KXYs(2,:,:));
KXs = sparse(KXs(:));
KYs = sparse(KYs(:));

parfor i = 1:numPoints*N
    % get components of beta
    bx = beta_xs(i);
    by = beta_ys(i);
    
    kxplusbx = KXs + bx;
    kyplusby = KYs + by;
    
    KX = diag(kxplusbx);
    KY = diag(kyplusby);

    % create e-value equation matrices
    AE = ( (KX/Ur)*KX + (KY/Ur)*KY  );
    AH = ( (KX/Er)*KX + (KY/Er)*KY  );
    
    BE = Er;
    BH = Ur;
    
     %Solve e-value equations
    %[VE,DE] = eig(AE,BE);
    %[VH,DH] = eig(AH,BH);
    
    AEsub_p = AE(SHInds_phys,SHInds_phys);
    BEsub_p = BE(SHInds_phys,SHInds_phys);
    AHsub_p = AH(SHInds_phys,SHInds_phys);
    BHsub_p = BH(SHInds_phys,SHInds_phys);
    
    AEsub_np = AE(SHInds_nonphys,SHInds_nonphys);
    BEsub_np = BE(SHInds_nonphys,SHInds_nonphys);
    AHsub_np = AH(SHInds_nonphys,SHInds_nonphys);
    BHsub_np = BH(SHInds_nonphys,SHInds_nonphys);

    % solve physical bands
    [VE_p,DE_p] = eig(AEsub_p,BEsub_p);
    [VH_p,DH_p] = eig(AHsub_p,BHsub_p);
    
    % Scale D
    DE_p = diag( sqrt(DE_p) * Lx / (2*pi) );
    DH_p = diag( sqrt(DH_p) * Lx / (2*pi) );
    
    % Normalise eigenvectors
    VE_p = normaliseEigVectors(VE_p);
    VH_p = normaliseEigVectors(VH_p);
    
    DE_dat_p(i,:) = DE_p;
    DH_dat_p(i,:) = DH_p;
    VE_dat_p(i,:,:) = VE_p;
    VH_dat_p(i,:,:) = VH_p;
    

    % solve non-physical bands
    [VE_np,DE_np] = eig(AEsub_np,BEsub_np);
    [VH_np,DH_np] = eig(AHsub_np,BHsub_np);
    
    % Scale D
    DE_np = diag( sqrt(DE_np) * Lx / (2*pi) );
    DH_np = diag( sqrt(DH_np) * Lx / (2*pi) );
    
    % Normalise eigenvectors
    VE_np = normaliseEigVectors(VE_np);
    VH_np = normaliseEigVectors(VH_np);
    
    DE_dat_np(i,:) = DE_np;
    DH_dat_np(i,:) = DH_np;
    VE_dat_np(i,:,:) = VE_np;
    VH_dat_np(i,:,:) = VH_np;
    
    AEs(i,:,:) = AE;
    AHs(i,:,:) = AH;
    BEs(i,:,:) = BE;
    BHs(i,:,:) = BH;
    
%     % Check which submatrix e-vectors belong to:
%     initialIndices = false(1,NH);
%     initialIndices(reshape((1:2:P)'+(0:1:Q-1)*P,[1,NSH])) = 1;
%     testMatE = abs(VE) < 1E-10;      % Find zero elements in E e-vectors
%     testMatH = abs(VH) < 1E-10;      % Find zero elements in H e-vectors
%     testVecE = (sum(testMatE(initialIndices,:)) + sum(~testMatE(~initialIndices,:))) < NSH; % Truth table for indices of primary submatrix e-vectors
%     testVecH = (sum(testMatH(initialIndices,:)) + sum(~testMatH(~initialIndices,:))) < NSH;
%     
%     DE1_dat(i,:) = DE(testVecE);
%     DH1_dat(i,:) = DH(testVecH);
end
fprintf([': ', num2str(toc(startBand)), 's\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorder eigenvalues and eigenvectors to same order for all betas %
DE_dat_sorted_p = sort(abs(DE_dat_p'))';     % Take the magnitude of the e-values - they should be real as matrix hermitian, but floating point errors give small imaginary component
DH_dat_sorted_p = sort(abs(DH_dat_p'))';

DE_dat_sorted_np = sort(abs(DE_dat_np'))';
DH_dat_sorted_np = sort(abs(DH_dat_np'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw Graphs %%

fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position = [5,5,18.5,18.5];
fig1.Color = [1 1 1];
pause(0.001);

bottom = 0;   % Find the highest and lowest value in the relative permittivities and permeabilities of both materials
top = max(max([eGrid,uGrid]));

% Plot Permittivity of unit cell
sp1 = subplot(3,3,1);
imagesc(eGrid);
title([char(949) '_r']);
caxis manual;
caxis([bottom top]);
sp1.XTick = [];
sp1.XTickLabel = [];
sp1.YTick = [];
sp1.YTickLabel = [];
colorbar;
pbaspect([Lx Ly 1]);
set(gca,'FontName', 'calibri');
set(gca,'FontSize', 12);

% Plot Permeability of unit cell
sp2 = subplot(3,3,4);
imagesc(uGrid);
title('\mu_r');
caxis manual;
caxis([bottom top]);
colorbar;
sp2.XTick = [];
sp2.XTickLabel = [];
sp2.YTick = [];
sp2.YTickLabel = [];
pbaspect([Lx Ly 1]);
set(gca,'FontName', 'calibri');
set(gca,'FontSize', 12);


% Plot betas explored in 1st BZ.
sp3 = subplot(3,3,7);
plot(b(1,:),b(2,:),'-');
hold on;
BZpoints = [ahash_+bhash_,ahash_-bhash_,-ahash_-bhash_,-ahash_+bhash_,ahash_+bhash_]/2;
plot(BZpoints(1,:),BZpoints(2,:),'k:');
title('\beta Trace');
pointPos = [points{:,1}];   %Get list of BZ point coordinates.
text(pointPos(1,:),pointPos(2,:),points(:,2))
maxrange = max([-pi/Lx, pi/Lx, -pi/Ly, pi/Ly]);
xlim([-maxrange,maxrange]*1.1);
ylim([-maxrange,maxrange]*1.1)
pbaspect([1 1 1]);
xticks([-pi,0,pi]);
xticklabels({'-\pi','0','\pi'});
yticks([-pi,0,pi]);
yticklabels({'-\pi ','0 ','\pi '});
set(gca,'FontName', 'calibri');
set(gca,'FontSize', 12);

% Plot the Band Structure
sp4 = subplot(3,3,[2,3,5,6,8,9]);
hold on

% Note, workaround for transparency in 2018 ( p.Color(4) = 0.5; ) no longer
% seems to work in 2023. Instead, using 'patchline'.
% To fix patches in legend, start with dummy plot
dummyplotH = plot([-200,-200],[-200,-200],'r-','LineWidth',2);
dummyplotE = plot([-200,-200],[-200,-200],'b-','LineWidth',2);

for t = (1:NSH_physical)
    %p = plot(dists(2:end),DH_dat_sorted_p(:,t),'r-'); % plot H/TE modes (red)
    p = patchline(dists(2:end),DH_dat_sorted_p(:,t),'edgecolor','r','linestyle','-','edgealpha',0.25); % plot H/TE modes (red)
    p.LineWidth = 2;
    %p.Color(4) = 0.5;
    %p = plot(dists(2:end),DE_dat_sorted_p(:,t),'b-'); % plot E/TM modes (blue)
    p = patchline(dists(2:end),DE_dat_sorted_p(:,t),'edgecolor','b','linestyle','-','edgealpha',0.25); % plot E/TM modes (blue)
    p.LineWidth = 2;
    %p.Color(4) = 0.5;
end
for t = (1:NSH_nonphysical)
    p = plot(dists(2:end),DH_dat_sorted_np(:,t),'r:'); % plot H/TE modes (red)
    %p = patchline(dists(2:end),DH_dat_sorted_np(:,t),'edgecolor','r','facecolor','r','linestyle',':','edgealpha',0.5); % plot H/TE modes (red)
    p.LineWidth = 2;
    %p.Color(4) = 0.5;
    p = plot(dists(2:end),DE_dat_sorted_np(:,t),'b:'); % plot E/TM modes (blue)
    %p = patchline(dists(2:end),DE_dat_sorted_np(:,t),'edgecolor','b','facecolor','b','linestyle',':','edgealpha',0.5); % plot E/TM modes (blue)
    p.LineWidth = 2;
    %p.Color(4) = 0.5;
end

axis([-0.01,dists(end)+0.01,0,0.8])
%title(['Band Structure of Infinite Square Crystal, ', num2str(P), ' Harmonics']);
xlabel('Bloch Wave Vector, \beta');
ylabel('Normalised Frequency');

% Get key point labels and dividing vertical lines plotted over band structure
xticks(dists((0:1:length(points))*N+1));
pointNamesFull = [points{:,2}];
pointNamesFull(length(points)+1) = points{1,2};
xticklabels(pointNamesFull);
line([1;1]*dists((0:1:length(points))*N+1) , [-1000;1000]*ones(length(points)+1,1)' , 'Color',[0 0 0],'LineStyle','--');

legend('H/TE modes','E/TM modes','Location','southeast');

set(gca,'FontName', 'calibri');
set(gca,'FontSize', 12);

hold off
fig = fig1;     % Set figure as output object.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Figure to file.
if ~isfolder('figures')
    mkdir('figures')
end

filename = 'figures//BS_';
filename = [filename,'_',char(datetime('now','Format','yyyy-MM-dd HHmm')),'.fig'];
saveas(fig,filename);

if closeFig
    close(fig);
end

endTime = tic();
fprintf(['Band Structure calculated in ', num2str(toc(startTime)), 's\n']);
