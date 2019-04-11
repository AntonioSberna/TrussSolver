function [elements] = elementsProperty(elements,nodes, sections, materials)
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
