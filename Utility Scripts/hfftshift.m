function X = hfftshift(x)
% hfft2 leaves the origin split across the ASA coordinate system as
%
%   2 1 0 0 0 0
%    1 0 0 0 0 0
%   0 0 0 0 0 1
%    0 0 0 0 0 1
%   0 0 0 0 0 0
%    1 0 0 0 0 0
%
% This function reorders the lowest frequency fourier components to be
% adjacent to eachother.

a0 = x.a0;
a1 = x.a1;

size0 = size(a0);
size1 = size(a1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 1 0 0 0 0 0 0           2 1 0 0 0 0 0 1 
%  1 0 0 0 0 0 0 0           1 0 0 0 0 0 0 1
% 0 0 0 0 0 0 0 0           0 0 0 0 0 0 0 0 
%  0 0 0 0 0 0 0 1  -----\   0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 1   -----/  0 0 0 0 0 0 0 0 
%  0 0 0 0 0 0 0 1           0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0           0 0 0 0 0 0 0 0 
%  1 0 0 0 0 0 0 0           1 0 0 0 0 0 0 1 

a0(:,floor(size0(2)/2)+1:size0(2)) = fftshift(a0(:,floor(size0(2)/2)+1:size0(2)),1);
a1(:,floor(size1(2)/2)+1:size1(2)) = fftshift(a1(:,floor(size1(2)/2)+1:size1(2)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 1 0 0 0 0 0 1 			0 0 0 0 0 0 0 0 
%  1 0 0 0 0 0 0 1			 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 			0 0 0 0 0 0 0 0 
%  0 0 0 0 0 0 0 0	-----\	 0 0 0 1 1 0 0 0 
% 0 0 0 0 0 0 0 0 	-----/	0 0 0 1 2 1 0 0 
%  0 0 0 0 0 0 0 0			 0 0 0 1 1 0 0 0 
% 0 0 0 0 0 0 0 0 			0 0 0 0 0 0 0 0 
%  1 0 0 0 0 0 0 1 			 0 0 0 0 0 0 0 0

a0 = fftshift(a0);
a1 = fftshift(a1);



X = HexImage(a0,a1);