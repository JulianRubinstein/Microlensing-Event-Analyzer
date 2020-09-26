clc; clear all; close all;

%input theoretical paramaters and input function
uo = 0.3;
to = 24;
te = 1;

%sampling the function at n different locations (boundaries here are "maximum area" that we want to fit the parabole to)
n = 50 ; %number of samples
a = to-0.03*te*(1-uo.^2) ; %left limit (dependent on te and uo, the relative width of the event)
b = to+0.03*te*(1-uo^2) ; %right limit (dependent on te uo, the relative width of the event)
t=linspace(a,b,n); %defining the vector space

%defining the theoritcal function
y = ( (uo).^2 + ((t-to)./(te)).^2 + 2 ) ./ ( sqrt((uo).^2 + ((t-to)./(te)).^2 ) .* sqrt((uo).^2 + ((t-to)./(te)).^2 + 4) );
dy = 0.005*ones(1,n); ; %defining an arbitrary error vector for simulation

%choosing number of samples for parabola paramaters and defining vectors
%for paramaters storage (parabola of the form (ax^2 + bx + c)).
%-a1 and b1 are not to be mixed with a and b wchich define the range for the
%parabola fit-

%para_to,para_uo,para_te are the approximations for the respective evenet
%paramaters which were approximated from the parabola.

sampnum=3000; %number of samples taken

c1=zeros(1,sampnum); % vector of c values
b1=zeros(1,sampnum); % vector of b values
a1=zeros(1,sampnum); % vector of a values
d_c1=zeros(1,sampnum); % vector of c error values
d_b1=zeros(1,sampnum); % vector of b error values
d_a1=zeros(1,sampnum); % vector of a error values
x_2=zeros(1,sampnum); % vector of x_2 values

%simulating data and fitting function sampnum times to extract paramaters
for j = 1:sampnum
    y1=noisesim(y,dy,n);
    [k,dk]=funcfit(y1,dy,t,n);
    
    %storing the values
    c1(j)=k(1);
    b1(j)=k(2);
    a1(j)=k(3);
    d_c1(j)=sqrt(dk(1,1));
    d_b1(j)=sqrt(dk(2,2));
    d_a1(j)=sqrt(dk(3,3));
    
    y_para=k(3)*t.^2+k(2).*t+k(1);
    x_2(j)=sum(((y1-y_para).^2)./((dy).^2));
end

%average value and standard deviation of a,b,c paramaters and errors
a1_avg=mean(a1);
b1_avg=mean(b1);
c1_avg=mean(c1);
a1_error=sqrt((std(a1)/sqrt(sampnum)).^2+(sqrt(1/(sum(1./d_a1.^2))))^2);
b1_error=sqrt((std(b1)/sqrt(sampnum)).^2+(sqrt(1/(sum(1./d_b1.^2))))^2); 
c1_error=sqrt((std(c1)/sqrt(sampnum)).^2+(sqrt(1/(sum(1./d_c1.^2))))^2);
x_2_avg=mean(x_2);
x_2_reduced=x_2_avg/(n-3);

%histograms for presenting data
figure('Name','a');
histogram(a1,floor(sqrt(sampnum)*2))
figure('name','b')
histogram(b1,floor(sqrt(sampnum)*2))
figure('name','c')
histogram(c1,floor(sqrt(sampnum)*2))

%analytic calculation of a,b,c with above paramaters for comparison
a_true=-17.91;
b_true=859.7112;
c_true=-10313.08961;