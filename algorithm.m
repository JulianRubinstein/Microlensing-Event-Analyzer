%clc; clear all; close all;

%This function accepts vectors y1 and t of length n which contain the magnitude
%of the effect (yi) for each time (xi). In addition it recieves initial guesses
%for uo, to and te. The following are such paramaters for testing:
% 
% uo = 0.5;
% to = 0;
% te = 2.5;
% n = 50 ;
% a = -0.3 ;
% b = 0.3 ;
% t=linspace(a,b,n);
% y = ( (uo).^2 + ((t-to)./(te)).^2 + 2 ) ./ ( ((uo).^2 + ((t-to)./(te)).^2 ) .* sqrt((uo).^2 + ((t-to)./(te)).^2 + 4) );
% dy = (0.01*ones(n,1))'; ; %defining an arbitrary error vector for simulation
% y1=zeros(1,n);
% for i=1:n
%    y1(i) = y(i) + normrnd(0,dy(i));
% end
% scatter(t,y1)
% 
% % Initial paramaters and defining theoretical function
% uo_exp=1.1357;
% to_exp=4.9986;
% te_exp=2.4969;
% d_uo=0.0037
% d_to=0.0167
% d_te=0.0037

%The actual function starts here:

function [final_x_2,final_uo,final_to,final_te,d_uo,d_to,d_te] = algorithm(y1,dy,t,uo_exp,to_exp,te_exp,d_uo,d_to,d_te)

init_uo=uo_exp;
init_to=to_exp;
init_te=te_exp;

y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) );
x_2=sum(((y1-y_exp).^2)./(dy.^2)); %initial chi square

%Defining chi squred and paramaters which are used for loops
sgn=1; %sign paramater which is used for flipping of the search direction
condition=0.001; %the accuracy in which we choose to find uo,to,te squared

%While loops which look for minimum of chi square in paramater space:
 
%first loop for "uo" paramater
epsilon=5*d_uo; %intial step size in paramter space
while(epsilon > condition)
    uo_exp = uo_exp+(sgn)*epsilon; %making a small step
    y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) ); %redefining the function after the small step
    x_2_new=sum(((y1-y_exp).^2)./(dy.^2)); %redefining chi square after small step
    if(x_2_new < x_2) %if chi square has been reduced update chi square
        x_2=x_2_new ;
    else %if chi square has increased, change the search direction and decrease the step size (as we are closer to minimum)
        sgn=(-1)*sgn;
        uo_exp = uo_exp+(sgn)*epsilon;
        epsilon = 0.5*epsilon;
    end
end

%same concept as loop #1 for "to" paramater
epsilon=5*d_to;
while(epsilon > 10*condition)
    to_exp = to_exp+(sgn)*epsilon;
    y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) );
    x_2_new=sum(((y1-y_exp).^2)./(dy.^2));
    if(x_2_new < x_2)
        x_2=x_2_new;
    else
        sgn=(-1)*sgn;
        to_exp = to_exp+(sgn)*epsilon;
        epsilon = 0.5*epsilon;
    end
end

%same concept as while #1 for "te" paramater
epsilon=5*d_te;
while(epsilon > condition)
    te_exp = te_exp+(sgn)*epsilon;
    y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) );
    x_2_new=sum(((y1-y_exp).^2)./(dy.^2));
    if(x_2_new < x_2)
        x_2=x_2_new;
    else
        sgn=(-1)*sgn;
        te_exp = te_exp+(sgn)*epsilon;
        epsilon = 0.5*epsilon;
    end
end

final_x_2=x_2;
final_uo=uo_exp;
final_te=te_exp;
final_to=to_exp;

%moving along the paramater space untill an increase of one chi_square has
%been detected, this is the uncertainity of the paramater.
epsilon=max(0.00001,(to_exp-init_to)/10);
while(abs(final_x_2-x_2)<1)
    to_exp=to_exp-epsilon;
    y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) );
    x_2=sum(((y1-y_exp).^2)./(dy.^2));
end
d_to=final_to-to_exp;
to_exp=final_to; %returning the paramater to the correct value
x_2=final_x_2; %returning chi square to the correct value

%same concept
epsilon=max(0.00001,(te_exp-init_te)/10);
while(abs(final_x_2-x_2)<1)
    te_exp=te_exp-epsilon;
    y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) );
    x_2=sum(((y1-y_exp).^2)./(dy.^2));
end
d_te=final_te-te_exp;
te_exp=final_te;
x_2=final_x_2;



epsilon=max(0.00001,(uo_exp-init_uo)/10);
while(abs(final_x_2-x_2)<1)
    uo_exp=uo_exp+epsilon;
    y_exp = ( (uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 2 ) ./ ( sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 ) .* sqrt((uo_exp).^2 + ((t-to_exp)./(te_exp)).^2 + 4) );
    x_2=sum(((y1-y_exp).^2)./(dy.^2));
end
d_uo=final_uo-uo_exp;
uo_exp=final_uo;
x_2=final_x_2;


