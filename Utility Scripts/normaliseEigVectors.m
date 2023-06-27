function Vout = normaliseEigVectors(Vin)
%% Takes a square matrix of column vectors, and returns the matrix of 
%  normalised column vectors, i.e. the modulus of each column vector is one

sizeVin = size(Vin);
onesCoVec = ones(sizeVin(1));

V2 = abs(Vin.^2);           % create matrix of the magnitude of each squared element.
V3 = sqrt(onesCoVec*V2);    % Create covector with magnitude of each eigVector.
Vout = Vin./V3;
