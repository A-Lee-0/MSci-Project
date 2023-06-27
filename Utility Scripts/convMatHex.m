function C = convMat(A,P,Q,R)
% CONVMAT Rectangular Convolution Matrix
%
% C = convmat(A,P); for 1D problems
% C = convmat(A,P,Q); for 2D problems
% C = convmat(A,P,Q,R); for 3D problems
%
% This MATLAB function constructs convolution matrices
% from a real-space grid.
%% HANDLE INPUT AND OUTPUT ARGUMENTS
% DETERMINE SIZE OF A
[N0x,N0y] = size(A.a0);
[N1x,N1y] = size(A.a1);
% HANDLE NUMBER OF HARMONICS FOR ALL DIMENSIONS
if nargin==2
Q = 1;
end

% COMPUTE INDICES OF SPATIAL HARMONICS
NH = P*Q; %total number
p = -floor(P/2):+floor(P/2); %indices along x
q = -floor(Q/2):+floor(Q/2); %indices along y


% COMPUTE FOURIER COEFFICIENTS OF A
%A = fftshift(fftn(A)) / (Nx*Ny*Nz);
A = hfftshift(hfft2(A));
A.a0 = A.a0./(N0x*N0y + N1x*N1y);
A.a1 = A.a1./(N0x*N0y + N1x*N1y);


% COMPUTE ARRAY INDICES OF CENTER HARMONIC
p0 = 1 + floor(N0x/2);
q0 = 1 + floor(N0y/2);


% INITIALIZE CONVOLUTION MATRIX
C = zeros(NH,NH);


% POPULATE CONVOLUTION MATRIX
for qrow = 1 : Q
for prow = 1 : P
    row = (qrow-1)*P + prow;
    for qcol = 1 : Q
    for pcol = 1 : P
        col = (qcol-1)*P + pcol;
        pfft = p(prow) - p(pcol);
        %disp(['pfft is equal to ',num2str(pfft),'.'])
        qfft = q(qrow) - q(qcol);
        %disp(['qfft is equal to ',num2str(qfft),'.'])
        
        %disp(['row is equal to ',num2str(row),'.'])
        %disp(['col is equal to ',num2str(col),'.'])
        
        % find hfft component
%         k1 = round(k_(1)*2/sqrt(3))/2;
%         k2 = round(k_(2)*2)/2;
%         if mod(k2,1) == 0
%             x = k1;
%             y = k2;
%             K = hexCellFFT.a0(501+x,251+y);
%         else
%             x = k1 - 0.5;
%             y = k2 - 0.5;
%             K = hexCellFFT.a1(501+x,251+y);
%         end
        
%         if mod(qfft,2) == 0
%             C(row,col) = A.a0(p0+pfft+(qfft/2),q0+ (qfft/2) );
%         else
%             C(row,col) = A.a1(p0+pfft +((qfft-1)/2),q0+ ((qfft-1)/2) );
%         end

        if mod(pfft,2) == 0
            C(row,col) = A.a0(p0-pfft/2, q0+qfft - pfft/2 );
        else
            C(row,col) = A.a1(p0-(pfft+1)/2, q0+qfft -(pfft+1)/2);
        end
        
    end
    end
end
end



