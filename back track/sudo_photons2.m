num_p = 10000000;
x = normrnd(0,0.6866,1,num_p);
y = normrnd(0,0.7615,1,num_p);

xb_po = [];
yb_po = [];

dmax = 14;
di = 0.05;
num_d = dmax/di;

pri_pho_i = [];
w = [];

d_dis = zeros(num_d,2);
for j = 1:num_d
    d_dis(j,1) = di*(j - 0.5); 
end

for i = 1:num_p    
    n = Ind_specific(i);
    x0 = x(i) + 26.7*par(n,4)/sqrt(1-(par(n,4)^2+par(n,5)^2));
    xb_po(end+1) = x0;
    y0 = y(i) + 26.7*par(n,5)/sqrt(1-(par(n,4)^2+par(n,5)^2));
    yb_po(end+1) = y0;
    d = sqrt(x0^2 + y0^2);

    for k = 1:num_d
        if di*(k-1) < d && d < di*k
            d_dis(k,2) = d_dis(k,2) + 1;
        end
    end
end

sum(d_dis(:,2))

figure(1)
x2 = d_dis(:,1);
y2 = d_dis(:,2)./d_dis(:,1);
scatter(x2,y2,'.');
xlabel('distance to origin/cm')
ylabel("relative number")
title('position distribution of back track photons on plane Z=26.7')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_p = 10000000;
x = normrnd(0,0.6866,1,num_p);
y = normrnd(0,0.7615,1,num_p);

xb_po = [];
yb_po = [];

dmax = 14;
di = 0.05;
num_d = dmax/di;

pri_pho_i = [];
w = [];

d_dis = zeros(num_d,2);
for j = 1:num_d
    d_dis(j,1) = di*(j - 0.5); 
end
for i = 1:num_p    
    n = Ind_specific(i);
    x0 = par(n,2) + 26.7*par(n,4)/sqrt(1-(par(n,4)^2+par(n,5)^2));
    xb_po(end+1) = x0;
    y0 = par(n,3) + 26.7*par(n,5)/sqrt(1-(par(n,4)^2+par(n,5)^2));
    yb_po(end+1) = y0;
    d = sqrt(x0^2 + y0^2);

    for k = 1:num_d
        if di*(k-1) < d && d < di*k
            d_dis(k,2) = d_dis(k,2) + 1;
        end
    end
end

figure(1)
x2 = d_dis(:,1);
y2 = d_dis(:,2)./d_dis(:,1);
scatter(x2,y2,'.');
xlabel('distance to origin/cm')
ylabel("relative number")
title('position distribution of back track photons on plane Z=26.7')
grid on