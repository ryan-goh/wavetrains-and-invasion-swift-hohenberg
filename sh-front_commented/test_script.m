%%test script to get continuation running;
clear all
close all

mu_init = 0.4;
gpuon = 0; 
pardirec = -[-1 1];



main_fcn(mu_init,gpuon,pardirec)