clc
clear all
close all
density=xlsread('my_results.xlsx',2);
velocity=xlsread('my_results.xlsx',1);
shear=xlsread('my_results.xlsx',3);
[~,idx] = unique(density(:,1));   %which rows have a unique first value?
density = density(idx,:);
[~,idx] = unique(velocity(:,1));   %which rows have a unique first value?
velocity = velocity(idx,:);
[~,idx] = unique(shear(:,1));   %which rows have a unique first value?
shear = shear(idx,:);
vden=interp1(density(:,1),density(:,2),shear(:,1));
vvel=interp1(velocity(:,1),velocity(:,2),shear(:,1));
f=zeros(38,1);
for i=1:38
f(i)=shear(i,2)/(1/2*vden(i)*vvel(i)*vvel(i));
end
plot (shear(:,1),f)
avef=mean(f)