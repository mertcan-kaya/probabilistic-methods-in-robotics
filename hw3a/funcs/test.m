clc, clear all, close all

n = 20000;
b_sq = 2;

nsum = zeros(1,23);

norm_dist = zeros(n,1);
for i = 1:n
    norm_dist(i,1) = sample_normal_distribution(b_sq);
    
    if norm_dist(i,1) < -5.25
        nsum(1) = nsum(1) + 1;
    elseif norm_dist(i,1) < -4.75
        nsum(2) = nsum(2) + 1;
    elseif norm_dist(i,1) < -4.25
        nsum(3) = nsum(3) + 1;
    elseif norm_dist(i,1) < -3.75
        nsum(4) = nsum(4) + 1;
    elseif norm_dist(i,1) < -3.25
        nsum(5) = nsum(5) + 1;
    elseif norm_dist(i,1) < -2.75
        nsum(6) = nsum(6) + 1;
    elseif norm_dist(i,1) < -2.25
        nsum(7) = nsum(7) + 1;
    elseif norm_dist(i,1) < -1.75
        nsum(8) = nsum(8) + 1;
    elseif norm_dist(i,1) < -1.25
        nsum(9) = nsum(9) + 1;
    elseif norm_dist(i,1) < -0.75
        nsum(10) = nsum(10) + 1;
    elseif norm_dist(i,1) < -0.25
        nsum(11) = nsum(11) + 1;
    elseif norm_dist(i,1) < 0.25
        nsum(12) = nsum(12) + 1;
    elseif norm_dist(i,1) < 0.75
        nsum(13) = nsum(13) + 1;
    elseif norm_dist(i,1) < 1.25
        nsum(14) = nsum(14) + 1;
    elseif norm_dist(i,1) < 1.75
        nsum(15) = nsum(15) + 1;
    elseif norm_dist(i,1) < 2.25
        nsum(16) = nsum(16) + 1;
    elseif norm_dist(i,1) < 2.75
        nsum(17) = nsum(17) + 1;
    elseif norm_dist(i,1) < 3.25
        nsum(18) = nsum(18) + 1;
    elseif norm_dist(i,1) < 3.75
        nsum(19) = nsum(19) + 1;
    elseif norm_dist(i,1) < 4.25
        nsum(20) = nsum(20) + 1;
    elseif norm_dist(i,1) < 4.75
        nsum(21) = nsum(21) + 1;
    elseif norm_dist(i,1) < 5.25
        nsum(22) = nsum(22) + 1;
    else
        nsum(23) = nsum(23) + 1;
    end
end

tsum = zeros(1,23);

trig_dist = zeros(n,1);
for i = 1:n
    trig_dist(i,1) = sample_triangular_distribution(b_sq);
    
    if trig_dist(i,1) < -5.25
        tsum(1) = tsum(1) + 1;
    elseif trig_dist(i,1) < -4.75
        tsum(2) = tsum(2) + 1;
    elseif trig_dist(i,1) < -4.25
        tsum(3) = tsum(3) + 1;
    elseif trig_dist(i,1) < -3.75
        tsum(4) = tsum(4) + 1;
    elseif trig_dist(i,1) < -3.25
        tsum(5) = tsum(5) + 1;
    elseif trig_dist(i,1) < -2.75
        tsum(6) = tsum(6) + 1;
    elseif trig_dist(i,1) < -2.25
        tsum(7) = tsum(7) + 1;
    elseif trig_dist(i,1) < -1.75
        tsum(8) = tsum(8) + 1;
    elseif trig_dist(i,1) < -1.25
        tsum(9) = tsum(9) + 1;
    elseif trig_dist(i,1) < -0.75
        tsum(10) = tsum(10) + 1;
    elseif trig_dist(i,1) < -0.25
        tsum(11) = tsum(11) + 1;
    elseif trig_dist(i,1) < 0.25
        tsum(12) = tsum(12) + 1;
    elseif trig_dist(i,1) < 0.75
        tsum(13) = tsum(13) + 1;
    elseif trig_dist(i,1) < 1.25
        tsum(14) = tsum(14) + 1;
    elseif trig_dist(i,1) < 1.75
        tsum(15) = tsum(15) + 1;
    elseif trig_dist(i,1) < 2.25
        tsum(16) = tsum(16) + 1;
    elseif trig_dist(i,1) < 2.75
        tsum(17) = tsum(17) + 1;
    elseif trig_dist(i,1) < 3.25
        tsum(18) = tsum(18) + 1;
    elseif trig_dist(i,1) < 3.75
        tsum(19) = tsum(19) + 1;
    elseif trig_dist(i,1) < 4.25
        tsum(20) = tsum(20) + 1;
    elseif trig_dist(i,1) < 4.75
        tsum(21) = tsum(21) + 1;
    elseif trig_dist(i,1) < 5.25
        tsum(22) = tsum(22) + 1;
    else
        tsum(23) = tsum(23) + 1;
    end
end

figure(1)
plot(norm_dist,0,'.')

figure(2)
bar(-5.5:0.5:5.5,nsum)

figure(3)
plot(trig_dist,0,'.')

figure(4)
bar(-5.5:0.5:5.5,tsum)
