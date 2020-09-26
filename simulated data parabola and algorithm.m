clc; clear all; close all;

%input theoretical paramaters and input function
uo = 0.3
to = 1
te = 1.5

%sampling the function at n different locations (boundaries here are "maximum area" that we want to fit the parabole to)
n = 50 ; %number of samples
a = to-0.2*te*(1-uo^2) ; %left limit (dependent on te and uo, the relative width of the event)
b = to+0.2*te*(1-uo^2) ; %right limit (dependent on te and uo, the relative width of the event)
t=linspace(a,b,n); %defining the vector space

%defining the theoritcal function
y = (( (uo).^2 + ((t-to)./(te)).^2 + 2 ) ./ ( sqrt((uo).^2 + ((t-to)./(te)).^2 ) .* sqrt((uo).^2 + ((t-to)./(te)).^2 + 4) ));
dy = 0.01*ones(1,n); %defining an arbitrary error vector for simulation

%choosing number of samples for parabola paramaters and defining vectors
%for paramaters storage (parabola of the form (ax^2 + bx + c)).
%a1 and b1 are not to be mixed with a and b wchich define the range for the
%parabola fit

sampnum=1000; %number of samples taken

te1=zeros(1,sampnum); % vector of te values
to1=zeros(1,sampnum); % vector of to values
uo1=zeros(1,sampnum); % vector of uo values
d_te1=zeros(1,sampnum); % vector of te error values
d_to1=zeros(1,sampnum); % vector of to error values
d_uo1=zeros(1,sampnum); % vector of uo error values
final_x_2_values=zeros(1,sampnum); %vector of x_2 values
final_x_2_reduced_values=zeros(1,sampnum); %vector of x_2 reduced values

%para_to,para_uo,para_te are the approximations for the respective evenet
%paramaters which were approximated from the parabola.

    y1=noisesim(y,dy,n);
    [k,dk]=funcfit(y1,dy,t,n);
    syms f(x)
    f(x)=k(3)*x^2+k(2)*x+k(1);
    
    fplot(f(x),[a,b])
    hold
    scatter(t,y1)

%simulating data and fitting function sampnum times to extract paramaters
for j = 1:sampnum
    
    y1=noisesim(y,dy,n);
    [k,dk]=funcfit(y1,dy,t,n);

    c1=k(1);
    b1=k(2);
    a1=k(3);
    c1_error=sqrt(dk(1,1));
    b1_error=sqrt(dk(2,2));
    a1_error=sqrt(dk(3,3));
    
    [para_to,para_uo,para_te,p,p_inv,d_experiment] = matrices(a1,b1,c1,a1_error,b1_error,c1_error);
    
    [final_x_2,final_uo,final_to,final_te,d_uo,d_to,d_te] = algorithm(y1,dy,t,para_uo,para_to,para_te,d_experiment(3),d_experiment(1),d_experiment(2));
    
    final_x_2_reduced=final_x_2/(n-3);
    
    %storing the values
    final_x_2_values(j)=final_x_2;
    final_x_2_reduced_values(j)=final_x_2_reduced;
    te1(j)=final_te;
    to1(j)=final_to;
    uo1(j)=final_uo;
    d_te1(j)=d_te;
    d_to1(j)=d_to;
    d_uo1(j)=d_uo;
end

%average value and errors of paramaters
x_2__avg=mean(final_x_2_values)
x_2_reduced_avg=mean(final_x_2_reduced_values)
te_avg=mean(te1)
to_avg=mean(to1)
uo_avg=mean(uo1)
te_error=sqrt((std(te1)/sqrt(sampnum)).^2+(sqrt(1/(sum(1./d_te1.^2))))^2)
to_error=sqrt((std(to1)/sqrt(sampnum)).^2+(sqrt(1/(sum(1./d_to1.^2))))^2) 
uo_error=sqrt((std(uo1)/sqrt(sampnum)).^2+(sqrt(1/(sum(1./d_uo1.^2))))^2) 
