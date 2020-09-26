function [k,dk] = funcfit(y,dy,t,n)

V=diag(dy.^2); %the covariance matrix (In this case, diagonal matrix of variances: V(i,i) = dy(i)^2)

%defining the matrix C of polynom powers (barlow 6.21) for parabole fit
c1=ones(n,1); %polynom of power 0 vector
c2=t'; %polynom of power 1 vector
c3=(t.^2)'; %polynom of power 2 vector
C = [c1, c2, c3]; %concatination

dk = ((C')*(V^(-1))*(C))^(-1);  %error matrix of paramaters; (barlow 6.24)
k = (dk)*((C')*V^(-1))*y'; %vector of length 3 which contains a,b,c (barlow 6.23)


