function [para_to,para_uo,para_te,p,p_inv,d_experiment] = matrices(a1,b1,c1,a1_error,b1_error,c1_error)

%parmaters g and d are for internal calclculations and are to be ignored

%defining the non-linear fucntion which transforms from a,b,c into uo,to,te.

g1=c1^2-(b1^2*c1/(2*a1))+b1^4/(16*a1^2);

para_to=-b1/(2*a1);
para_uo=abs(sqrt((2*(((1-g1)+sqrt(g1^2-g1))./(g1-1)))));
para_te=abs((8*para_to/(b1*para_uo^3*(para_uo^2+4)^(3/2)))^(1/2));

%defining p matrix which transforms from d_uo,d_to_d_te to d_a, d_b, d_c
%and p_inv which does the opposite.

d1=-8/(para_uo^3*(para_uo^2+4)^(3/2)*para_te^2);
d2=16/(para_uo^3*(para_uo^2+4)^(3/2)*para_te^3);
d3=(48*para_uo^2+96)/(para_te^2*para_uo^4*(para_uo^2+4)^(5/2));

p = zeros(3);
p(1,1) = 0;
p(1,2) = (d2/2).^2;
p(1,3) = (d3/2).^2;
p(2,1) = (-d1).^2;
p(2,2) = (-d2*para_to).^2;
p(2,3) = (d3*para_to).^2;
p(3,1) = (d1*para_to).^2;
p(3,2) = (-d2*((para_to^2)/2)).^2;
p(3,3) = (d3*((para_to^2)/2)-8/(para_uo^2*(para_uo^2+4)^(3/2))).^2 ;
p_inv = p^(-1);

%defining vector ((d_a^2,d_b^2,d_c^2) to transform in to (d_to^2,d_te^2,d_uo^2))
d_para=[a1_error.^(2),b1_error.^(2),c1_error.^(2)]';
d_experiment=sqrt(abs(p_inv*d_para));