num_p = 100000000;
x = normrnd(0,0.6866,1,num_p);
y = normrnd(0,0.7615,1,num_p);
% beam_divergence = normrnd(0,0.001,1,num_p);
% % w = cos(beam_divergence./2);
% w = cos(beam_divergence);
% theta = rand()*2*pi;
% m = sqrt(1-w.^2);
% u = m*cos(theta);
% v = m*cos(theta);
% 
% for i = 1:100
%     u(i)^2+v(i)^2+w(i)^2
% end



x1 = [];
y1 = [];
for i = 1:1000
    x1(end+1) = x(i) + u(i)/w(i)*26.7;
    y1(end+1) = y(i) + v(i)/w(i)*26.7;
end

scatter(x1,y1)

dmax = 1;
dmin = -1;
di = 0.005;
num_d = (dmax-dmin)/di;

d_dis = zeros(num_d,2);
for j = 1:num_d
    d_dis(j,1) = di*(j-0.5)+dmin; 
end

for i = 1:length(y1)
    x0 = y1(i);
    for k = 1:length(d_dis(:,1))
        if di*(k-1)+dmin < x0 && x0 < di*k+dmin
            d_dis(k,2) = d_dis(k,2) + 1;
        end
    end
end

figure()
% x2 = xb_po;
% y2 = yb_po;

x2 = d_dis(:,1);
y2 = d_dis(:,2);

xq = 0:0.000001:2;
vq = interp1(x2,y2,xq);
[~,in_half] = min(abs(vq-0.5*max(y2)));
FWHM = xq(in_half)

plot(x2,y2,'.');
xlabel('distance to origin/cm')
ylabel("relative number")
title('position distribution of back track photons on plane Z=0')
% xlabel('x/cm')
% ylabel("y/cm")
% title('position distribution of back track photons on plane Z=0')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xb_po = [];
yb_po = [];
dmax = 1;
dmin = -1;
di = 0.005;
num_d = (dmax-dmin)/di;

d_dis = zeros(num_d,2);
for j = 1:num_d
    d_dis(j,1) = di*(j-0.5)+dmin; 
end

for i = 1:num_what(what)     
    n = Ind_specific(i);
    x0 = par(n,2) - 26.7*par(n,4)/sqrt(1-(par(n,4)^2+par(n,5)^2));
    xb_po(end+1) = x0;
    y0 = par(n,3) - 26.7*par(n,5)/sqrt(1-(par(n,4)^2+par(n,5)^2));
    yb_po(end+1) = y0;

    for k = 1:length(d_dis(:,1))
        if di*(k-1)+dmin < y0 && y0 < di*k+dmin
            d_dis(k,2) = d_dis(k,2) + 1;
        end
    end
end

figure()
% x2 = xb_po;
% y2 = yb_po;

x3 = d_dis(:,1);
y3 = d_dis(:,2);

xq = 0:0.000001:2;
vq = interp1(x2,y2,xq);
[~,in_half] = min(abs(vq-0.5*max(y2)));
FWHM = xq(in_half)

plot(x2,y2,'.');
xlabel('distance to origin/cm')
ylabel("relative number")
title('position distribution of back track photons on plane Z=0')
% xlabel('x/cm')
% ylabel("y/cm")
% title('position distribution of back track photons on plane Z=0')
grid on