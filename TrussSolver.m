
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                        %%%   
%%%                  Automatic calculus                    %%%
%%%                     Homework n.4                       %%%
%%%                                                        %%%
%%%                       Group 16                         %%%
%%%                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all
clear 
clc

tic;
%% Part 1 - Input data

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

% Length and inclination of each frame

for i=1:size(elements,1)
    % Coordinates of the nodes of the i^th element
    x = [nodes(find(nodes(:,1)==elements(i,2),1),2) nodes(find(nodes(:,1)==elements(i,3),1),2)];
    y = [nodes(find(nodes(:,1)==elements(i,2),1),3) nodes(find(nodes(:,1)==elements(i,3),1),3)];  
    
    % Length of i^th element
    elements(i,5) = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);
    % Tilt angle of i^th elements
    elements(i,6) = atan((y(2) - y(1))/(x(2) - x(1)));
    % Area of i^th element
    elements(i,7) = sections(find(sections(:,1) == elements(i,4),1),8);
    % Elasticity modulus of the i^th element (element 4 -> section 9 -> material 2)
    elements(i,8) = materials(find(materials(:,1) == sections(find(sections(:,1) == elements(i,4),1),9),1),2);    
end
clear i j k x_1 x_2 y_1 y_2

%% Original structure plot break

orig_structural_plotter(nodes, restraints, forces, elements)
uiwait



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



% Deformed shape magnification factor
mag = 10;





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

%% Part 5 - Plot

structural_plotter(nodes, d, restraints, forces, elements, sigma, disp)


%% Part 6 - Outputs export

% Nodes output
Node(1,1) = "ID";
x_undef(1,1) = "[m]" ;
y_undef(1,1) = "[m]";
x_def(1,1) = "[m]";
y_def(1,1) = "[m]";

for i = 1:size(nodes,1)
    Node(i+1,1) = num2str(nodes(i,1));
    x_undef(i+1,1) = num2str(nodes(i,2));
    y_undef(i+1,1) = num2str(nodes(i,3));
    x_def(i+1,1) = num2str(round(def(i,2),2));
    y_def(i+1,1) = num2str(round(def(i,3),2));
end
NODES_OUTPUT = table(Node,x_undef, y_undef, x_def, y_def);

writetable(NODES_OUTPUT,'RESULTS_01_Nodes.txt','Delimiter','\t')
clear Node x_undef y_undef x_def y_def NODES_OUTPUT i


% Elements output
Element(1,1) = "ID";
Stress(1,1) = "[MPa]";
Strain(1,1) = "[%]";
f_ix(1,1) = "[kN]";
f_iy(1,1) = "[kN]";
f_jx(1,1) = "[kN]";
f_jy(1,1) = "[kN]";


for i = 1:size(nodes,1)
    Element(i+1,1) = num2str(nodes(i,1));
    Stress(i+1,1) = num2str(round(sigma(i,2),2));
    Strain(i+1,1) = num2str(round(epsilon(i,2),4)*100);
    f_ix(i+1,1) = num2str(f{i}((elements(i,2)*2-1),1)/1000);
    f_iy(i+1,1) = num2str(f{i}((elements(i,2)*2),1)/1000);
    f_jx(i+1,1) = num2str(f{i}((elements(i,3)*2-1),1)/1000);
    f_jy(i+1,1) = num2str(f{i}((elements(i,3)*2),1)/1000);
end
ELEMENTS_OUTPUT = table(Element, Stress, Strain, f_ix, f_iy, f_jx, f_jy);

writetable(ELEMENTS_OUTPUT,'RESULTS_02_Elements.txt','Delimiter','\t')
clear Element sigmas epsilons f_ix f_iy f_jx f_jy ELEMENTS_OUTPUT i


toc;



%% Attachment 1 - Plot functions

function structural_plotter(nodes, d, restraints, forces, elements, sigma, disp)
    %%%
    %%% commento da rivedere 
    %%%
    %%% Plot a structural system composed of "nodes", "restraints" (if res = 1)
    %%% labels near the node (if "label" = 1 and there is not a restrain), 
    %%% forces in kN (if "forc" = 1) and beams (according to "conn" connectivity 
    %%% matrix. It also uses magnification factors (separatelly for x and y directions.
    %%% As "colour" it needs a rgb vector.
    %%%
    %%% This funcion needs "restraints_draw" and "forces_drawer" functions.
    %%%


%% First parameters
% Colours of the structures
orig_col= [0.8 0.8 0.8];
def_col = [0 0 0];

% Stress colours legend
stresscolour = jet;


% First case magnification factor
mag = 10;
% Nodal coordinates of deformated shape
def_shape = zeros(size(nodes,1),2);
var = 1;
for i = 1:2:size(d)
    def_shape(var,1) = var;
    def_shape(var,2) = d(i)*mag + nodes(var,2);
    def_shape(var,3) = d(i+1)*mag + nodes(var,3);
    var = var + 1;
end
clear i var



%% Figure 
% Full screen windows
fig = figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1],'Visible','off');

hold on

% Magnification factors
x_max = max([max(def_shape(:,2)) max(nodes(:,2))]);
y_max = max([max(def_shape(:,3)) max(nodes(:,3))]);
x_min = min([min(def_shape(:,2)) min(nodes(:,2))]);
y_min = min([max(def_shape(:,3)) min(nodes(:,3))]);

magn_fact_x = (x_max+x_min)/20;
magn_fact_y = (y_max+y_min)/10;


%% Option panel

	p = uipanel(fig,'Title','Plot options','Position',[.1 .1 .1 .25]);


% Save button
    save_push = uicontrol(fig,'Style','pushbutton',...
        'String','Save Image','Callback',@savepushbutton_Callback,'Units','normalized', 'Position', [.21 .1 .08 .04]);
    saving = uicontrol(fig,'Style','text',...
        'String','Saving...','Units','normalized', 'Position', [.2 .05 .1 .04],'Visible','off');


% Check boxes
        elem_label = uicontrol(p,'Style','checkbox','Value',0,...
            'String','Elements n.','Callback',@elemlabel_Callback,'Units','normalized', 'Position', [.2 .8 .8 .15]);
        node_label = uicontrol(p,'Style','checkbox','Value',1,...
            'String','Nodes n.','Callback',@nodelabel_Callback,'Units','normalized', 'Position', [.2 .6 .8 .15]);
        displval = uicontrol(p,'Style','checkbox','Value',0,...
            'String','Displacements','Callback',@displval_Callback,'Units','normalized', 'Position', [.25 .45 .8 .15]);
        stresscol = uicontrol(p,'Style','checkbox','Value',1,...
            'String','Stress colour','Callback',@stresscol_Callback,'Units','normalized', 'Position', [.2 .25 .8 .15]);
        stressval = uicontrol(p,'Style','checkbox','Value',0,...
            'String','Stress values','Callback',@stressval_Callback,'Units','normalized', 'Position', [.25 .1 .8 .15]);
        




%% Colour Legend
%Thick labels for colour legend
thicklab = zeros(11,1);
for i = 1:11
thicklab(i) = round(min(sigma(:,2)) + (i-1)*((max(sigma(:,2)) - min(sigma(:,2)))/10),0);
end

%Colour legend
colormap(stresscolour);
col_bar = colorbar('TickLabels',{ num2str(thicklab(1)), num2str(thicklab(2)), num2str(thicklab(3)), ...
    num2str(thicklab(4)), num2str(thicklab(5)), num2str(thicklab(6)), num2str(thicklab(7)), ... 
    num2str(thicklab(8)), num2str(thicklab(9)), num2str(thicklab(10)), num2str(thicklab(11))},...
    'Location', 'eastoutside');
col_bar_MPa = text(1.08,0,sprintf(' MPa'), 'HandleVisibility','off', 'FontSize', 8,'Units','normalized');


%% Non deformed shape
% Non-deformed beams
for i = 1:size(elements,1)
    eval('orig_beam_%i',i) = plot([nodes(elements(i,2),2) nodes(elements(i,3),2)],...
        [nodes(elements(i,2),3) nodes(elements(i,3),3)], 'Color', orig_col);
end
clear i

% Non-deformed nodes
for i = 1:size(nodes,1)
    eval('orig_nodes_%i',i) = scatter(nodes(i,2),nodes(i,3),[], orig_col, 'HandleVisibility','off','MarkerFaceColor',[1 1 1]);
end
clear i




%% Deformed shape
% Restraints
for i = 1:size(restraints,1)
        % Hinge
    if restraints(i,2) == 1 && restraints(i,3) == 1
        restrain_type = 2;
        % Roller
    elseif restraints(i,2) == 1 || restraints(i,3) == 1
         restrain_type = 1;
    end
    restraints_drawer(def_shape(restraints(i,1),2),def_shape(restraints(i,1),3),...
        magn_fact_x,magn_fact_y,restrain_type,def_col);
end
clear i



% Deformed shape beams
for i = 1:size(elements,1)
    eval('def_beam_%i',i) = plot([def_shape(elements(i,2),2) def_shape(elements(i,3),2)],...
        [def_shape(elements(i,2),3) def_shape(elements(i,3),3)], 'Color', rainbow_plot(i,sigma(:,2),stresscolour), 'LineWidth',2);
end
clear i
  

% Deformed nodes
for i = 1:size(def_shape,1)
    eval('def_nodes_%i',i) = scatter(def_shape(i,2),def_shape(i,3),[], def_col, 'HandleVisibility','off','MarkerFaceColor',[1 1 1]);
end
clear i
    

% Forces
magn = (magn_fact_x+magn_fact_y)/2;
for i=1:size(forces,1)
    forces_drawer(def_shape, forces(i,:), magn, def_col)
end
clear i


% Magnification factor
str = sprintf('Magnification factor of deformed shape = %i',mag);
title(str)
set(get(gca,'title'),'Position',[(x_max+x_min)/2 y_min-6*magn_fact_y], 'FontSize',8)



%% Labels
% Labels near the elements
x = zeros(size(elements,1),1);
y = zeros(size(elements,1),1);
for i = 1:size(elements,1)
    x(i) = (def_shape(find(def_shape(:,1) == elements(i,2),1),2) + def_shape(find(def_shape(:,1) == elements(i,3),1),2))/2 + magn_fact_x*sin(elements(i,5));
    y(i) = (def_shape(find(def_shape(:,1) == elements(i,2),1),3) + def_shape(find(def_shape(:,1) == elements(i,3),1),3))/2 + magn_fact_y*cos(elements(i,5))/2;
end
t_elem_labels = text(x,y,num2str(elements(:,1)), 'HandleVisibility', 'off','Visible','off');



% Labels tension near the elements
x_label = zeros(size(elements,1),1);
y_label = zeros(size(elements,1),1);
for i = 1:size(elements,1)
    x_label(i) = (def_shape(find(def_shape(:,1) == elements(i,2),1),2) + def_shape(find(def_shape(:,1) == elements(i,3),1),2))/2 - magn_fact_x*sin(elements(i,5))/2;
    y_label(i) = (def_shape(find(def_shape(:,1) == elements(i,2),1),3) + def_shape(find(def_shape(:,1) == elements(i,3),1),3))/2 - magn_fact_y*cos(elements(i,5))/2;
end
t_sigma_labels = text(x_label,y_label, num2str(round(sigma(:,2))) + " MPa", 'HandleVisibility', 'off','Visible','off');



% Labels near the nodes 
node_labels = num2str(def_shape(:,1));
t_node_labels = text(def_shape(:,2)+magn_fact_x/2,def_shape(:,3)+ magn_fact_y/2,node_labels, 'HandleVisibility', 'off');


%% Plot dimentions
% ricontrollare
ylim([y_min-10*magn_fact_y y_max+5*magn_fact_y]) 
xlim([x_min-4*magn_fact_x x_max+4*magn_fact_x])

% Other parameters
axis off
fig.Visible = 'on';





%% Figure options

% Checkbox element labels
function elemlabel_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        set(t_elem_labels,'Visible','on');
    else
        set(t_elem_labels, 'Visible','off');
    end
end


% Checkbox node labels
function nodelabel_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        set(t_node_labels,'Visible','on')
        set(displval,'Visible','on')
    else
        set(t_node_labels,'Visible','off');
        set(displval,'Visible','off')
    end
end


% Checkbox stress colours
function stresscol_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value') == get(hObject,'Max'))
        % Deleting old beams
        for i = 1:size(elements,1)
            delete(eval('def_beam_%i',i));
        end
        clear i
        % Coloured deformed shape beams according tension state
        for  i = 1:size(elements,1)
           eval('def_beam_%i',i) = plot([def_shape(elements(i,2),2) def_shape(elements(i,3),2)],...
               [def_shape(elements(i,2),3) def_shape(elements(i,3),3)], 'Color', rainbow_plot(i,sigma(:,2),stresscolour), 'LineWidth',2);
        end
        clear i 
        
        set(col_bar,'Visible','on')
        set(col_bar_MPa,'Visible','on')
        set(stressval,'Visible','on')
        
    else
        % Deleting old beams
        for i = 1:size(elements,1)
            delete(eval('def_beam_%i',i));
        end
        clear i
        
        % Non - tinted deformed shape beams
        for  i = 1:size(elements,1)
           eval('def_beam_%i',i) = plot([def_shape(elements(i,2),2) def_shape(elements(i,3),2)],...
               [def_shape(elements(i,2),3) def_shape(elements(i,3),3)], 'Color', def_col); 
        end
        clear i 
        set(col_bar,'Visible','off')
        set(col_bar_MPa,'Visible','off')
        set(stressval,'Visible','off')
    end
end

% Checkbox stress values
function stressval_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        set(t_sigma_labels,'Visible','on');
    else
        set(t_sigma_labels,'Visible','off');
    end
end

% Checkbox diplacements values
function displval_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value') == get(hObject,'Max'))
        delete(t_node_labels)
        node_labels = num2str(def_shape(:,1)) + " ["  + num2str(round(disp(:,1),2)) +"," + num2str(round(disp(:,2),2)) +  "] m";
        t_node_labels = text(def_shape(:,2)+magn_fact_x/2,def_shape(:,3)+ magn_fact_y/2,node_labels, 'HandleVisibility', 'off');
    else
        delete(t_node_labels)
        node_labels = num2str(def_shape(:,1));
        t_node_labels = text(def_shape(:,2)+magn_fact_x/2,def_shape(:,3)+ magn_fact_y/2,node_labels, 'HandleVisibility', 'off');
    end
    
end

% Push button to save the images
push_counter = 1;
function savepushbutton_Callback(hObject,eventdata,handles)
    % Preparation to save
    set(p,'Visible','off')
    set(saving,'Visible','on')
    set(saving,'String',sprintf('Saving TrussSolver_%i.png',push_counter))
    
    % Save as png
    print(fig,sprintf('TrussSolver_%i',push_counter),'-dpng','-noui')

    % After save option
    pause(0.5)
    push_counter = push_counter + 1;
    set(saving,'String','Done!')
    pause(0.3)
    set(saving,'Visible','off')
    set(p,'Visible','on')

end



end

function orig_structural_plotter(nodes, restraints, forces, elements)
    %%%
    %%% commento da rivedere 
    %%%
    %%% Plot a structural system composed of "nodes", "restraints" (if res = 1)
    %%% labels near the node (if "label" = 1 and there is not a restrain), 
    %%% forces in kN (if "forc" = 1) and beams (according to "conn" connectivity 
    %%% matrix. It also uses magnification factors (separatelly for x and y directions.
    %%% As "colour" it needs a rgb vector.
    %%%
    %%% This funcion needs "restraints_draw" and "forces_drawer" functions.
    %%%


%% First parameters
% Colours of the structures
orig_col= [0 0 0];


% Magnification factors
x_max = max(nodes(:,2));
y_max = max(nodes(:,3));
x_min = min(nodes(:,2));
y_min = min(nodes(:,3));

magn_fact_x = (x_max+x_min)/20;
magn_fact_y = (y_max+y_min)/10;

%% Figure 
% Full screen windows
fig = figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1],'Visible','off');

hold on

%% Option panel

	p = uipanel(fig,'Title','Plot options','Position',[.1 .1 .1 .15]);


% Save button
    ok_push = uicontrol(fig,'Style','pushbutton',...
        'String','Proceed','Callback',@doitpushbutton_Callback,'Units','normalized', 'Position', [.21 .1 .08 .04]);
    calculating = uicontrol(fig,'Style','text',...
        'String',"Let's do it!",'Units','normalized', 'Position', [.2 .05 .1 .04],'Visible','off');


% Check boxes
        elem_label = uicontrol(p,'Style','checkbox','Value',1,...
            'String','Elements n.','Callback',@elemlabel_Callback,'Units','normalized', 'Position', [.2 .7 .8 .25]);
        node_label = uicontrol(p,'Style','checkbox','Value',1,...
            'String','Nodes n.','Callback',@nodelabel_Callback,'Units','normalized', 'Position', [.2 .25 .8 .25]);

%% Non deformed shape
% Non-deformed beams
for i = 1:size(elements,1)
    eval('orig_beam_%i',i) = plot([nodes(elements(i,2),2) nodes(elements(i,3),2)],...
        [nodes(elements(i,2),3) nodes(elements(i,3),3)], 'Color', orig_col);
end
clear i

% Restraints
for i = 1:size(restraints,1)
        % Hinge
    if restraints(i,2) == 1 && restraints(i,3) == 1
        restrain_type = 2;
        % Roller
    elseif restraints(i,2) == 1 || restraints(i,3) == 1
         restrain_type = 1;
    end
    restraints_drawer(nodes(restraints(i,1),2),nodes(restraints(i,1),3),...
        magn_fact_x,magn_fact_y,restrain_type,orig_col);
end
clear i


% Non-deformed nodes
for i = 1:size(nodes,1)
    eval('orig_nodes_%i',i) = scatter(nodes(i,2),nodes(i,3),[], orig_col,...
        'HandleVisibility','off','MarkerFaceColor',[1 1 1]);
end
clear i



    

% Forces
magn = (magn_fact_x+magn_fact_y)/2;
for i=1:size(forces,1)
    forces_drawer(nodes, forces(i,:), magn, orig_col)
end
clear i



%% Labels
% Labels near the elements
x = zeros(size(elements,1),1);
y = zeros(size(elements,1),1);
for i = 1:size(elements,1)
    x(i) = (nodes(find(nodes(:,1) == elements(i,2),1),2) + nodes(find(nodes(:,1) == elements(i,3),1),2))/2 + magn_fact_x*sin(elements(i,5));
    y(i) = (nodes(find(nodes(:,1) == elements(i,2),1),3) + nodes(find(nodes(:,1) == elements(i,3),1),3))/2 + magn_fact_y*cos(elements(i,5))/2;
end
t_elem_labels = text(x,y,num2str(elements(:,1)), 'HandleVisibility', 'off','Visible','on');



% Labels near the nodes 
node_labels = num2str(nodes(:,1));
t_node_labels = text(nodes(:,2)+magn_fact_x/2,nodes(:,3)+ magn_fact_y/2,node_labels, 'HandleVisibility', 'off');




% Magnification factor
str = sprintf("This is not the output, click 'Proceed' to run the code!");
title(str)
set(get(gca,'title'),'Position',[(x_max+x_min)/2 y_min-7*magn_fact_y], 'FontSize',11,'Color',[1 .1 0])


%% Plot dimentions
% ricontrollare
ylim([y_min-10*magn_fact_y y_max+5*magn_fact_y]) 
xlim([x_min-4*magn_fact_x x_max+4*magn_fact_x])

% Other parameters
axis off
fig.Visible = 'on';


%% Figure options

% Checkbox element labels
function elemlabel_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        set(t_elem_labels,'Visible','on');
    else
        set(t_elem_labels, 'Visible','off');
    end
end


% Checkbox node labels
function nodelabel_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        set(t_node_labels,'Visible','on')
    else
        set(t_node_labels,'Visible','off');
    end
end


% Push button to save the images
push_counter = 1;
function doitpushbutton_Callback(hObject,eventdata,handles)
    % Preparation to save
    set(p,'Visible','off')
    set(calculating,'Visible','on')
    pause(0.6)
    close(fig)
end
end 
    
function [] = restraints_drawer(x,y,fact_x,fact_y,type,colour)
    %%%
    %%% 1 for rollers
    %%% 2 for hinges
    %%% joints not jet implemented
    %%%

    % Points coordinates
    if type == 1    % roller case
        a = [x-fact_x y-(fact_y*1.5)];
        b = [x-(fact_x/2) y-fact_y];
        c = [x y];
        d = [x+(fact_x/2) y-fact_y];
        e = [x+fact_x y-(fact_y*1.5)];
    elseif type == 2    % hinge case
        a = [x-fact_x y-fact_y];
        b = [x-(fact_x/2) y-fact_y];
        c = [x y];
        d = [x+(fact_x/2) y-fact_y];
        e = [x+fact_x y-fact_y];
    elseif type == 3 % horizontal joint case
    elseif type == 4 % vertical joint case
    end
    
    
    % Lines coordinates
    xs = [b(1) c(1) d(1) b(1) a(1) e(1)];
    ys = [b(2) c(2) d(2) b(2) a(2) e(2)];
    
    
    % Point plot
    plot(xs(1:4),ys(1:4), 'Color', colour, 'HandleVisibility', 'off')
    plot(xs(5:6),ys(5:6), 'Color', colour, 'HandleVisibility', 'off')
    
    %Field plot
    y_field = [a(2) a(2)-(fact_y/2)];
    field_lenght = (xs(6) - xs(5))/4;
    for i = 0:4
        x_field = [xs(5)+field_lenght*i xs(5)+field_lenght*(i-1)];
        plot(x_field, y_field,'Color', colour, 'HandleVisibility', 'off')
    end 
    clear i
    

    % Roller wheels
    if type == 1
        r = (b(2)-a(2))/2;
        y = b(2) - r;
        x1 =(b(1) + (d(1) - b(1))/4);
        x2 =(d(1) - (d(1) - b(1))/4);
        circle(x1,y,r,colour)
        circle(x2,y,r,colour)
    end

    function circle(x,y,r,colour)
        circ=0:0.01:2*pi; 
        xc=r*cos(circ);
        yc=r*sin(circ);
        plot(x+xc,y+yc, 'Color', colour, 'HandleVisibility','off');
            
    end
end


function [] = forces_drawer(node, force, magn, colour)
    %%%
    %%% Plot force (attenzione alla forma dei vettori di input!)
    %%%
    
    %% Horizontal forces
    i = 2;
    if force(i) ~= 0
        p1 = [(node(force(1),2)-(force(i)/(2000000*magn))) node(force(1),3)];
        p2 = [node(force(1),2) node(force(1),3)];
        dp = p2-p1;
        quiver(p1(1),p1(2),dp(1),dp(2),0, 'Color', colour)
        text((p1(1)+p2(1)-magn)/2,(p1(2)+p2(2)+magn)/2, sprintf('%.0f kN',force(i)/1000)) 
    end
    clear i    
    
    %% Vertical forces
    i = 3;
    if force(i) ~= 0
        p1 = [node(force(1),2) (node(force(1),3)-(force(i)/(2000000*magn)))];
        p2 = [node(force(1),2) node(force(1),3)];
        dp = p2-p1;
        quiver(p1(1),p1(2),dp(1),dp(2),0, 'Color', colour)
        text((p1(1)+p2(1))/2,(p1(2)+p2(2)+magn)/2, sprintf('%.0f kN',abs(force(i))/1000)) 
    end 
    clear i
  %https://it.mathworks.com/matlabcentral/answers/160487-how-can-i-draw-a-line-with-arrow-head-between-2-data-points-in-a-plot
end

% Rainbow scale
function col = rainbow_plot(i, values,colormap)
    index = round(interp1([min(values),max(values)],[1,64],values(i)));
    col = colormap(index,:);
end  
