%fhplot plots free hand test data and draws a circle encompassing 60% of the data points.
%d = fhplot(a, tilt, [center])
%a = bivarate scatter data
%tilt = tilt angle in degrees
%encomp = percentage of your data that you would like encompassed by a
%           circle
%lim = angular limit for plot
%center = center of circle drawn around 60% of the data. If blank, center will be set to the centroid.

function d = fhplot(a, tilt, encomp, lim,center)

x = (a(:,3)-lim);
y = (a(:,4)-lim);
[m,n] = size(x);
[m,n] = size(y);
f = sum(x)/m; %finds the center of mass for the polar coord scatter
p = sum(y)/m;

if nargin < 4, center = [f,p]; end

for i = 1:m;
	distance(i) = sqrt((x(i)-center(1))^2 + (y(i)-center(2))^2); %finding the distance between all the pionts and the center
end

dist = sort([distance]); %sorting the distances into an array
[c e] = size(dist); 
num = round((encomp/100)*e);
r = dist(round(num)) %setting the radius of your circle to a percentage of the total number of points

tilt2 = tilt - (tilt/2);
tilt3 = tilt + (tilt/2);
tilt4 = tilt3 + (tilt/2);
figure1 = figure;
axes1 = axes('Parent',figure1,'DataAspectRatio',[lim lim 1]);
box(axes1,'on')
plot(x,y,'k*','Parent',axes1)
hold on
ylim(axes1,[-lim lim])
ylim(axes1,[-lim lim])
plot(center(1),center(2),'ro','LineStyle','none','Color',[1 0 0])
theta = linspace(0, 2*pi); %generates a fxn "theta" that forms a complete circle
plot(r*cos(theta)+(center(1)), r*sin(theta)+(center(2)),'r','LineWidth',2,'Color',[1 0 0]) %circle denoting your angular error
plot(tilt*cos(theta), tilt*sin(theta),'Color',[0 0 0]) %circle denoting your applied tilt
plot(tilt2*cos(theta), tilt2*sin(theta),'Color',[0 0 0]);
plot(tilt3*cos(theta), tilt3*sin(theta),'Color',[0 0 0]);
plot(tilt4*cos(theta), tilt4*sin(theta),'Color',[0 0 0]);
yaxis = linspace(-lim,lim);
xaxis = linspace(-lim,lim);
plot(0,yaxis);
plot(xaxis,0);
xlabel('Tilt direction (degrees)','FontWeight','bold','FontSize',14,...
    'FontName','Times New Roman');
ylabel('Tilt direction (degrees)','FontWeight','bold','FontSize',14,...
    'FontName','Times New Roman');
str = sprintf('%02d percent of particles are within %02d degrees of expected angle',encomp,round(r));
title({str},...
    'FontSize',14,...
    'FontName','Times');
hold off
end