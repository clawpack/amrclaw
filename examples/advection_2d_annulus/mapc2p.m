function [xp,yp] = mapc2p(xc,yc)
% Map to polar coordinates.

xp = xc.*cos(yc);
yp = xc.*sin(yc);
