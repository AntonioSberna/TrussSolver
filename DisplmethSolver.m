function [nodes, d, restraints, forces, elements, sigma, disp, epsilon] = DisplmethSolver(materials, sections, nodes, elements, restraints, forces)

%% Part 3 - Stiffness matrix

K = zeros(size(nodes,1)*2);
k_loc_all = cell(1,size(elements,2)-1);
k_glob_all = cell(1,size(elements,2)-1);

for i=1:size(elements,1)
    % Sine and cosine of element tilt angle
    m=sin(elements(i,6));
    l=cos(elements(i,6));
    
    %Direction cosine matrix of the i^th element
    T = [l*l m*l -l*l -m*l;m*l m*m -m*l -m*m;-l*l -l*m l*l l*m;-l*m -m*m l*m m*m];
    % Local stiffness matrix
    k_loc = (elements(i,8)*elements(i,7)/elements(i,5))*T;
    k_loc_all{i} = k_loc;
    
    % Assemblage matrix
    A = zeros(size(nodes,1)*2,size(k_loc,1));
    A(2*elements(i,2)-1,1) = 1;
    A(2*elements(i,2),2) = 1;
    A(2*elements(i,3)-1,3) = 1;
    A(2*elements(i,3),4) = 1;
        
    % Assemblage
    k_glob = A*k_loc*A';
    k_glob_all{i} = k_glob;
    K = K + k_glob;
end
clear i k_loc k_glob m l A T

%% Part 4 - Displacement

% Nodal forces vector

F = zeros(size(nodes,1)*2,1);

for i = 1:size(forces,1)
  F(forces(i,1)*2-1,1) =  forces(i,2);
  F(forces(i,1)*2,1) =  forces(i,3);
end
clear i



% Reduction of the stiffness matrix (in order to figure out the
% unknown displacement values)

K_red = K;
F_red = F;

d_cont = ones(size(nodes,1)*2,1);
for i = size(restraints,1): -1 :1
  for j = 3:-1:2
    if restraints(i,j) == 1
      if j == 3
        K_red(restraints(i,1)*2,:) = [];
        K_red(:,restraints(i,1)*2) = [];
        F_red(restraints(i,1)*2) = [] ;
        d_cont(restraints(i,1)*2) = 0;
      else
        K_red(restraints(i,1)*2-1,:) = [];
        K_red(:,restraints(i,1)*2-1) = [];
        F_red(restraints(i,1)*2-1) = [];
        d_cont(restraints(i,1)*2-1) = 0;
      end
    end
  end
end
clear i j 


% Ammissible displacement
d_red = K_red\F_red;


% Expansion of the displacement vector
var = 1;
d = zeros(size(nodes,1)*2,1);
for i = 1:size(nodes,1)*2
  if d_cont(i) == 1
    d(i) = d_red(var);
    var = var+1;
  end
end
clear i var d_cont F_red K_red d_red



% Calculation of complete force vector (including the reaction values)
F = K * d;


% Nodal forces 
f = cell(1,size(elements,2)-1);
for i = 1:length(k_glob_all)
    f{i} = k_glob_all{i} * d;
end
clear i


% Nodal coordinates of deformated shape
def = zeros(size(nodes,1),2);
var = 1;
for i = 1:2:size(d)
    def(var,1) = var;
    def(var,2) = d(i) + nodes(var,2);
    def(var,3) = d(i+1) + nodes(var,3);
    var = var + 1;
end
clear i var



% Rearrangement of the displacement vector
disp = zeros(size(nodes,1),2);
var = 1;
for i = 1:2:size(d)
    disp(var,1) = d(i);
    disp(var,2) = d(i+1);
    var = var + 1;
end
clear i var


% Strains and stresses inside the bars
sigma = zeros(size(elements,1),2);
epsilon = zeros(size(elements,1),2);
L_def = zeros(size(elements,1),1);

for i = 1:size(elements,1)   
%   d_temp = [displ(find(elements(i,2) == displ(:,1),1),2)
%                displ(find(elements(i,2) == displ(:,1),1),3)
%                displ(find(elements(i,3) == displ(:,1),1),2)
%                displ(find(elements(i,3) == displ(:,1),1),3)];

%   c = cos(elements(i,6));
%   s = sin(elements(i,6));


%   sigma(i) = (elements(i,8) / elements(i,5))* [-c -s c s] * d_temp * (10^-6);
%   epsilon(i) = sigma(i) / elements(i,8);


   L_def(i) = ((def(find(nodes(:,1) == elements(i,2),1),2) - def(find(nodes(:,1) == elements(i,3),1),2))^2 ...
                 + (def(find(nodes(:,1) == elements(i,2),1),3) - def(find(nodes(:,1) == elements(i,3),1),3))^2)^0.5;
   
   epsilon(i,1) = elements(i,1);
   sigma(i,1) = elements(i,1);
   
   epsilon(i,2) = (L_def(i)-elements(i,5))/(elements(i,5));
   sigma(i,2) = elements(i,8) * epsilon(i,2)*(10^-6);
   
   


   
end
clear i c s


 %% Part 4 - Output export

txtexport(nodes,def,sigma,epsilon,f,elements)