function structural_plotter(nodes, d, restraints, forces, elements, sigma, disp, epsilon, def, f)
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
    saveim_push = uicontrol(fig,'Style','pushbutton',...
        'String','Save Image','Callback',@saveimpushbutton_Callback,'Units','normalized', 'Position', [.21 .1 .08 .04]);
    saving = uicontrol(fig,'Style','text',...
        'String','Saving...','Units','normalized', 'Position', [.25 .05 .1 .04],'Visible','off');
    savetxt_push = uicontrol(fig,'Style','pushbutton',...
        'String','Save Image','Callback',@savetxtpushbutton_Callback,'Units','normalized', 'Position', [.3 .1 .08 .04]);

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
function saveimpushbutton_Callback(hObject,eventdata,handles)
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


% Push button to save the images
function savetxtpushbutton_Callback(hObject,eventdata,handles)
    % Preparation to save
    set(p,'Visible','off')
    set(saving,'Visible','on')
    set(saving,'String','Saving results in txt file')
    
    % Save txt results file
    txtexport(nodes,def,sigma,epsilon,f,elements)

    % After save option
    pause(0.5)
    set(saving,'String','Done!')
    pause(0.3)
    set(saving,'Visible','off')
    set(p,'Visible','on')

end



end


function col = rainbow_plot(i, values,colormap)

% Rainbow scale

    index = round(interp1([min(values),max(values)],[1,64],values(i)));
    col = colormap(index,:);
end  
