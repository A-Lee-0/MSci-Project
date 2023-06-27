function b = DrawHexImage(x)
% DrawHexImage takes a HexImage structure, and puts the data into a single
% composite array to be used with the imagesc function. returns the
% composite array.

size0 = size(x.a0);
nx = size0(2);
ny = size0(1);

b = zeros(2*ny,2*nx+1);
for i = 1:nx
    for j = 1:ny
        b(2*j - 1, 2*i - 1) = x.a0(j,i);
        b(2*j - 1, 2*i) = x.a0(j,i);
        b(2*j, 2*i) = x.a1(j,i);
        b(2*j, 2*i + 1) = x.a1(j,i);
    end
end

