%%Appendix B

% Draw the planet
function DrawEarth(R)
[xx, yy, zz] = sphere(100);
surf(R*xx, R*yy, R*zz);
r = 0.8; g = r; b = r;
map = [r g b
       0 0 0
       r g b];
colormap(map)
caxis([-R/100 R/100])
shading interp
% Draw and label the X, Y and Z axes
line([0 2*R], [0 0], [0 0]); text(2*R, 0, 0, 'X')
line( [0 0], [0 2*R], [0 0]); text( 0, 2*R, 0, 'Y')
line( [0 0], [0 0], [0 2*R]); text( 0, 0, 2*R, 'Z')