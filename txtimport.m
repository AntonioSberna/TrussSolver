function [materials, sections, nodes, elements, restraints, forces] = txtimport()

% Import materials
file = fopen('01_Materials.txt','r');
temp = textscan(file, '%f %f %f %f' , 'headerLines' , 2);

materials = zeros(length(temp{1,1}),length(temp));
for i=1:length(temp)
    materials(:,i) = cell2mat(temp(i));  % ( ID , E , G , alpha )
end
clear i temp file

% Import sections
file = fopen('02_Sections.txt','r');
temp = textscan(file, '%s %f %f %f %f %f %f %f %f %f' , 'headerLines' , 2);

% N.B. POTREMO PENSARE DI TOGLIERE LA PRIMA COLONNA CON IL NUMERO DI
% SEZIONE, NON SERVE
sections = zeros(length(temp{1,1}),length(temp)-1);
for i=2:1:length(temp)
    sections(:,i-1) = cell2mat(temp(i));    
end
clear i temp file


% Import coordinates
file = fopen('03_Coordinates.txt','r');

temp = textscan(file, '%f %f %f' , 'headerLines' , 2);

nodes = zeros(length(temp{1,1}),length(temp));
for i=1:length(temp)
    nodes(:,i) = cell2mat(temp(i)); % ( ID , x , y )
end
clear i temp file


% Import connectivity
file = fopen('04_Elements.txt','r');

temp = textscan(file, '%f %f %f %f' , 'headerLines' , 1);

elements = zeros(length(temp{1,1}),length(temp));
for i=1:length(temp)
    elements(:,i) = cell2mat(temp(i)); % ( Element , Node_i , Node_j , length , ID_Section )
end
clear i temp file


% Import restraint
file = fopen('05_Restraints.txt' , 'r');

temp = textscan(file, '%f %f %f' , 'headerLines' , 2);

restraints = zeros(length(temp{1,1}),length(temp));
for i=1:length(temp)
    restraints(:,i) = cell2mat(temp(i)); % ( ID , T_x , T_y )
end
clear i temp file


% Import external forces
file = fopen('06_Forces.txt' , 'r');

temp = textscan(file, '%f %f %f' , 'headerLines' , 2);

forces = zeros(length(temp{1,1}),length(temp));
for i=1:length(temp)
    forces(:,i) = cell2mat(temp(i)); % ( ID , F_x , F_y )
end
clear i temp file




%% Part 2 - Charateristic of each element

[elements] = elementsProperty(elements,nodes, sections, materials);


%% Part 3 - Original structure plot break

orig_structural_plotter(nodes, restraints, forces, elements)
uiwait
end
