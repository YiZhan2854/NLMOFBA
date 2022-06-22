%code for base case 
%%%%%%%%%%%%%%%%%%%%%%%%%set parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_one=1; 
%ATP production rate for glycolysis
m_three=15; 
%ATP production rate for respiration
alpha_one=0.5; 
% enzyme cost for glycolysis 
alpha_two=0.5; 
% enzyme cost for fermentation 
alpha_three=10; 
% enzyme cost for respiration 
totalcost=200; 
% total enzyme resource
V_glucose=10000; 
% glucose availability upper bound
X_one=10000; 
% reaction rate upper bound for glycolysis
X_three=10000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_one=2;
c_two=1;
c_three=7; %normalization factors 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matrix to record results for plotting
percentageOfFermentation=zeros(52,52); 
resp=zeros(52,52); 
ferm=zeros(52,52); 
glyco=zeros(52,52);
ATPrate=zeros(52,52);
yield=zeros(52,52);
lactate=zeros(52,52);
%multi-objective optimization 
j=1;
for a=0:0.02:1
    k=1;
    for b=0:0.02:1
        c=1-a-b;
        if c<-eps %negative c is not realistic, set these values for plot
            percentageOfFermentation(j,k)=-1; 
            resp(j,k)=-10; 
            ferm(j,k)=-100; 
            glyco(j,k)=-100;
            ATPrate(j,k)=-200;
            yield(j,k)=-10;
            lactate(j,k)=-100;
            k=k+1;
            continue
        end    
%discretize the feasible region
        v_one_right=totalcost;
        v_one_left=0;
        v_one_interval=zeros(1,2002);
        for i=1:1:2001
            v_one_interval(i)=0.1*(i-1);
        end
        v_one_interval(2002)=totalcost/(alpha_one+alpha_three);
        v_two_lowerbound=zeros(1,2002);
        optimalvalues=zeros(1,4004);
        v_two_upperbound=v_one_interval;
%calculate values on the feasible region
        for i=1:1:2002
            optimalvalues(i)=c_one*a*m_one*v_one_interval(i)+c_two*b*v_two_upperbound(i)+c_one*a*m_three*(v_one_interval(i)-v_two_upperbound(i))-2*(1-a-b)*c_three*m_three*(v_two_upperbound(i)/v_one_interval(i));
        end
        for i=1:1:2002
            if v_one_interval(i)<totalcost/(alpha_one+alpha_three)
                v_two_lowerbound(i)=0;
            else
                v_two_lowerbound(i)=((alpha_one+alpha_three)/(alpha_three-alpha_two))*v_one_interval(i)-totalcost/(alpha_three-alpha_two);
            end
            optimalvalues(i+2002)=c_one*a*m_one*v_one_interval(i)+c_two*b*v_two_lowerbound(i)+c_one*a*m_three*(v_one_interval(i)-v_two_lowerbound(i))-2*(1-a-b)*c_three*m_three*(v_two_lowerbound(i)/v_one_interval(i));
        end
% find the optimal
        tempcopy=optimalvalues;
        [optimal,Index] = max(tempcopy);
        tempcopy(Index)=tempcopy(Index)-10000; 
        [optimal_sec,Index_sec] = max(tempcopy);
% determine the uniqueness of the optimal point
        if (abs(optimal-optimal_sec)<10*eps) && (~(v_one_interval(mod(Index,2002))==0.1 && v_one_interval(mod(Index_sec,2002))==0.2)) && (~(v_one_interval(mod(Index,2002))==totalcost && v_one_interval(mod(Index_sec,2002))==totalcost))  
            disp('more than one optimal point!');
        end
        if Index<2003
            glyco(j,k)=v_one_interval(Index);
            ferm(j,k)=v_two_upperbound(Index);
% continue determining the uniqueness of the optimal point
            if abs(c_two*b-c_one*a*m_three-2*(1-a-b)*c_three*m_three/v_one_interval(Index))<10*eps
                fprintf('optimal points could be along one entire line!');
            end
        else
            glyco(j,k)=v_one_interval(Index-2002);
            ferm(j,k)=v_two_lowerbound(Index-2002);
% continue determining the uniqueness of the optimal point
            if abs(c_two*b-c_one*a*m_three-2*(1-a-b)*c_three*m_three/v_one_interval(Index-2002))<10*eps
                fprintf('optimal points could be along one entire line!');
            end
        end
% compute some values of interest 
        resp(j,k)=glyco(j,k)-ferm(j,k); 
        percentageOfFermentation(j,k)=ferm(j,k)/glyco(j,k); 
        ATPrate(j,k)=m_one*glyco(j,k)+m_three*resp(j,k);
        yield(j,k)=2*ATPrate(j,k)/glyco(j,k);
        lactate(j,k)=ferm(j,k);
        k=k+1;
    end
    j=j+1;
end
    

%figure 1 ATP rate
figure(1);
pcolor(ATPrate)
color = gray;
color = flipud(color);
colormap(color);
colorbar;
title("ATP Production Rate");
ylabel("Weight of ATP production rate"); 
xlabel("Weight of lactate generation"); 
xticks([1 6 11 16 21 26 31 36 41 46 51]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
yticks(1:5:51);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
%figure 2 fermentation/(respiration+fermentation) extent of WE
figure(2);
pcolor(percentageOfFermentation)
color = gray;
color = flipud(color);
colormap(color);
colorbar;
title("fermentation/(respiration+fermentation)");
ylabel("Weight of ATP production rate"); 
xlabel("Weight of lactate generation");
xticks([1 6 11 16 21 26 31 36 41 46 51]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
yticks(1:5:51);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
%figure 3 fermentation
figure(3);
pcolor(ferm)
color = gray;
color = flipud(color);
colormap(color);
colorbar;
title("fermentation rate or lactate production rate");
ylabel("Weight of ATP production rate"); 
xlabel("Weight of lactate generation");
xticks([1 6 11 16 21 26 31 36 41 46 51]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
yticks(1:5:51);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
%figure 4 respiration
figure(4);
pcolor(resp)
color = gray;
color = flipud(color);
colormap(color);
colorbar;
title("respiration rate");
ylabel("Weight of ATP production rate"); 
xlabel("Weight of lactate generation");
xticks([1 6 11 16 21 26 31 36 41 46 51]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
yticks(1:5:51);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
%figure 5 glycolysis
figure(5);
pcolor(glyco)
color = gray;
color = flipud(color);
colormap(color);
colorbar;
title("glycolysis rate");
ylabel("Weight of ATP production rate"); 
xlabel("Weight of lactate generation");
xticks([1 6 11 16 21 26 31 36 41 46 51]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
yticks(1:5:51);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
%figure 6 yield
figure(6);
pcolor(yield)
color = gray;
color = flipud(color);
colormap(color);
colorbar;
title("yield");
ylabel("Weight of ATP production rate"); 
xlabel("Weight of lactate generation");
xticks([1 6 11 16 21 26 31 36 41 46 51]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
yticks(1:5:51);
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
