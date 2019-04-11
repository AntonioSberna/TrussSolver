function [] = txtexport(nodes,def,sigma,epsilon,f,elements)
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