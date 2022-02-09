
% File contain two function

%    1) for objective function (decide target condition with the help compound variable )
%    2) for developing system of first order equations using compound matrix method.

function target_bc_value=rec_plate_objective_func(l)

xmesh=linspace(-1,1);                  % size of rectangular plate (integration limit) 

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

[a b]=ode45(@(x,w)rec_plate_compound(x,w,l),xmesh,[0 1 0 0 0 0],options); %%%% Intial condition corresponds to boundary condition W'(-1)=W'''(-1)=0

target_bc_value=(b(end,5))^2;   % target condition achieved taking determinant of (CM) at right boundary 

end



%% This function converts fourth order system (W'''' + psi_2 W'' + psi_0 W ) into sixth order compound matrix system

function dwdx=rec_plate_compound(x,w,l)
global h 
% n=7;
alpha1=0.2;
psi_0=alpha1/(2*h);
psi_2=(l^4-1)*(2+(6+h*alpha1)*l^4)/(l^2+3*l^6);
psi_4=(4*h^2)/(3+9*l^4)*(3+h*alpha1+(2+3*h*alpha1)*l^4+(3+2*h*alpha1)*l^8);

%%% These equation are coming from paper (doi: 10.1016/j.jmps.2017.10.017)
A=[0                1    0              0;
   0                0    1              0;
   0                0    0              1;
   -psi_0/psi_4     0   -psi_2/psi_4    0];

%%% Following equations are derived using compound matrix method.

        f1=w(2);
  f2=w(3)+w(4);
  f3=A(4,2)*w(1)+A(4,3)*w(2)+A(4,4)*w(3)+w(5);
  f4=w(5);
  f5=-A(4,1)*w(1)+A(4,3)*w(4)+A(4,4)*w(5)+w(6);
  f6=-A(4,1)*w(2)-A(4,2)*w(4)+A(4,4)*w(6);
  
  dwdx=[f1; f2; f3; f4; f5; f6]; 
  
end
  


