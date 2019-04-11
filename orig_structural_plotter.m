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