clear; clc; close all;
load('FULL_SIM+ERR_hw_data_Thu_Mar_11_14-41-02_2021.mat')
t = timestamp(10:80);
% % plot(t,u,'b',t,v,'r')
states = x_hist(:,10:80);
clearvars -except t states
u = states(7,:);
v = states(8,:);
plot(t,u,'b',t,v,'r')
