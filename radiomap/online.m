close all, clear all, clc;


global M; %number of RSS beacons
global K; %number of TOA beacons

global R; %communication range of TOA beacons
global Delta; %variance of TOA measurement
global Npoint;


global Area; %sensing field;
global RssBeaconCoordinates; %coordinates of RSS beacons
global TOABeaconCoordinates; %coordinates of TOA beacons

global Gridsize;
global Map;
global Map_beacon;
global T; %training time (long enough)
global T_run; %running time
global H; %RSS comparing threshold

load beacon_setup.mat;

Y=zeros(M,M);
T_run = 10;
H=0;

col=Area(1)/GridSize;
row=Area(2)/GridSize;

Y_update = Y_offline;

%training location coefficients, size: 10*96000
%compare begin&update with accurate to demonstrate update is better than
%begin, so that to prove its adaptive property
X_accurate = []; %accurate coefficients, calculated in each time slot
X_begin = [];    %coefficients at time t=1
X_update = [];   %coefficients if updated

target = [];

for t=1:T_run %system running  
    %calculate current time's RSS by PropModel
    Rss = PropModel(Dis);
    Rss_location = PropModel(Dis_location);
    %generate target motion curve
    x = Area(1)/T_run*t;
    target(t,1)=x;
    target(t,2)=Area(2)/2*(sin(x*pi*2/Area(1))+1);
    %compute vector RSS_target, target<->beacon
    target_temp = [target(t,1),target(t,2)];
    Rss_target = Compute_Rss(target_temp, RssBeaconCoordinates);
    
    for i=1:M
        for j=1:M
            if i==j
                %do not calculate, ignore
            else
                %reconstruct Aij and Qij
                Q_ij = [];
                for k = 1:t
                    Q_ij = [Q_ij; Rss(i,j)];
                end
                A_ij = [];
                for k = 1:t
                    rss=[];
                    for k = 1:M
                        if k==i || k==j
                            
                        else
                            rss = [rss, Rss(i,k)];
                        end
                    end
                    rss = [rss, 1];
                    A_ij = [A_ij; rss];
                end
                temp = M*(i-1)+j;
                Y_ij = Y_offline(temp,:)';
                Q_temp = A_ij * Y_ij;
                %compare calculate RSS to accurate RSS
                %show whether to update Y_ij
                abs = (Q_temp - Q_ij).^2;             
                if sqrt(sum(sum(abs))) > H
                    Y_ij = inv(A_ij'*A_ij)*A_ij'*Q_ij; 
                    %replace in Y_update
                    Y_update(temp,:) = Y_ij;
                    %use Y_update to update training locations's coefficients
                    
                end
            end
        end
    end
    %calculate location coefficients first 
    %
    if t == 1
        X_begin = Cal_coefficient(Rss_location,Rss);
        X_update = X_begin;
    else
        %use t=1 X_begin to update training locations's coefficients adaptively
        %X_update = Cal_update_coefficient(Rss_location,Rss,X_begin,Y_update,M);
        %Update X_update at time t-1, generate New X_uodate at time t
        X_update = Cal_update_coefficient(Rss_location,Rss,X_update,Y_update,M);
    end
    %
    %calculate each time slot's accurate coefficients
    X_accurate = Cal_coefficient(Rss_location,Rss);
    %compare X_begin, update and accurate
    %want to show update is better than begin for its adaptive property
    flag_compare = 0;
    flag_compare = Compare_accuracy(X_begin,X_update,X_accurate);
    
    % simulate a motion curve target to verify
    flag_target = 0;
    flag_target = Compute_location(X_begin,X_update,target_temp,Rss_target,Rss);
end













