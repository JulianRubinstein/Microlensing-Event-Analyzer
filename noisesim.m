function [y1] = noisesim(y,dy,n)

y1=zeros(1,n);

for i=1:n
    y1(i) = y(i) + normrnd(0,dy(i));
end


