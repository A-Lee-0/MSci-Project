% Rumpf method adapted to use hexagonal basis, and hexagonal fft fourier
% coefficients.

% Calculates the band structure for a 2D hexagonal photonic crystals using
% the hexagonal FFT.
% It operates on the ASA coordinate structure eGrid and uGrid, containg the
% crystal unit cell's εᵣ and μᵣ respectively.

% This script calculates the band structure for a simple dielectric
% material with circular air holes drilled in a regular, hexagonal pattern.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup

startTime = tic();
pool = gcp();

% Define lattice vectors
a1_ = [1; 0; 0];
a2_ = [0.5; sqrt(3)/2; 0];

Lx = norm(a1_);
Ly = norm(a2_);

P = 13; 
Q = 13;
NH = P*Q; % total number of harmonics

N = 200;        % The number of bloch vectors to calculate between points in reciprocal space.

er = 13;        % The permittivity of the dielectric material
e0 = 1;         % The permittivity of the air
radius = 0.48;  % The radius of the air holes, as a fraction of the hexagon in-circle

resX = 500;     % The resolution of the unit cell, in X and Y. Note that in ASA, both offset arrays take this resolution!
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

% Construct Unit Cell.
a0 = ones(resY,resX)*er;    % populate default ASA arrays.
a1 = ones(resY,resX)*er;
rd = radius*resY;

c0_ = [-resX * sqrt(3) /2; 0];      % Centres of circles to define unit cell
c1_ = [-resX * sqrt(3)/2; resY];
c2_ = [resX * sqrt(3)*3/2; 0];
c3_ = [resX * sqrt(3)*3/2; resY];
c4_ = [resX * sqrt(3)/2; resY/2];

% loop through ASA arrays, and set any coordinates within circles as air.
for i = 1:resX
    for j = 1:resY
        r_ = [sqrt(3)*i; j];
        if norm(r_ - c0_) < rd
            a0(j,i) = e0;
        elseif norm(r_ - c1_) < rd
            a0(j,i) = e0;
        elseif norm(r_ - c2_) < rd
            a0(j,i) = e0;
        elseif norm(r_ - c3_) < rd
            a0(j,i) = e0;
        elseif norm(r_ - c4_) < rd
            a0(j,i) = e0;
        end
        
        r_ = r_ + [sqrt(3)/2; 1/2];
        if norm(r_ - c0_) < rd
            a1(j,i) = e0;
        elseif norm(r_ - c1_) < rd
            a1(j,i) = e0;
        elseif norm(r_ - c2_) < rd
            a1(j,i) = e0;
        elseif norm(r_ - c3_) < rd
            a1(j,i) = e0;
        elseif norm(r_ - c4_) < rd
            a1(j,i) = e0;
        end
    end
end

eGrid = HexImage(a0,a1);
uGrid = HexImage(ones(resY,resX),ones(resY,resX));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute convolution matrices %%
fprintf('Computing Convolution Matrices');
startConv = tic();

Er = convMatHex(eGrid,P,Q);
Ur = convMatHex(uGrid,P,Q);

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
points = [{[0;0],"\Gamma"}; 
          {bhash_/2,"M"}; 
          {(ahash_ + 2 * bhash_)/3,"K"}];
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

VE_dat = zeros(numPoints*N,P*Q,P*Q);        % Arrays to store the eigenvectors 
VH_dat = zeros(numPoints*N,P*Q,P*Q);
DE_dat = zeros(numPoints*N,P*Q);
DH_dat = zeros(numPoints*N,P*Q);

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
    [VE,DE] = eig(AE,BE);
    [VH,DH] = eig(AH,BH);
    
    % Scale D
    DE = diag( sqrt(DE) * Lx / (2*pi) );
    DH = diag( sqrt(DH) * Lx / (2*pi) );
    
    % Normalise eigenvectors
    VE = normaliseEigVectors(VE);
    VH = normaliseEigVectors(VH);
    
    DE_dat(i,:) = DE;
    DH_dat(i,:) = DH;
    VE_dat(i,:,:) = VE;
    VH_dat(i,:,:) = VH;
end
fprintf([': ', num2str(toc(startBand)), 's\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorder eigenvalues and eigenvectors to same order for all betas %
DE_dat_sorted = sort(abs(DE_dat'))';     % Take the magnitude of the e-values - they should be real as matrix hermitian, but floating point errors give small imaginary component
DH_dat_sorted = sort(abs(DH_dat'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw Graphs %%

fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position = [5,5,18.5,18.5];
fig1.Color = [1 1 1];
pause(0.001);

bottom = 0;   % Find the highest and lowest value in the relative permittivities and permeabilities of both materials
top = max(max([eGrid.a0,eGrid.a1,uGrid.a0,uGrid.a1]));

% Plot Permittivity of unit cell
sp1 = subplot(3,3,1);
imagesc(DrawHexImage(eGrid));
title([char(949) '_r']);
caxis manual;
caxis([bottom top]);
sp1.XTick = [];
sp1.XTickLabel = [];
sp1.YTick = [];
sp1.YTickLabel = [];
colorbar;
pbaspect([resX*sqrt(3) resY 1]);
set(gca,'FontName', 'calibri');
set(gca,'FontSize', 12);

% Plot Permeability of unit cell
sp2 = subplot(3,3,4);
imagesc(DrawHexImage(uGrid));
title('\mu_r');
caxis manual;
caxis([bottom top]);
colorbar;
sp2.XTick = [];
sp2.XTickLabel = [];
sp2.YTick = [];
sp2.YTickLabel = [];
pbaspect([resX*sqrt(3) resY 1]);
set(gca,'FontName', 'calibri');
set(gca,'FontSize', 12);

% Plot betas explored in 1st BZ.
sp3 = subplot(3,3,7);
plot(b(1,:),b(2,:),'-');
hold on;
BZpoints = [ahash_+2*bhash_,2*ahash_+bhash_,ahash_-bhash_,-ahash_-2*bhash_,-2*ahash_-bhash_,-ahash_+bhash_,ahash_+2*bhash_]/3;
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
for t = (1:NH)
    p = plot(dists(2:end),DH_dat_sorted(:,t),'r-'); % plot H/TE modes (red)
    p.LineWidth = 2;
    p = plot(dists(2:end),DE_dat_sorted(:,t),'b-'); % plot E/TM modes (blue)
    p.LineWidth = 2;
end
axis([-0.01,dists(end)+0.01,0,0.8])
%title(['Band Structure of Infinite Hexagonal Crystal, ', num2str(P), ' Harmonics']);
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
