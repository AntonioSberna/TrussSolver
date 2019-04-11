
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                        %%%   
%%%                      Truss Solver                      %%%
%%%                                                        %%%
%%%                                                        %%%
%%%            antoniopio.sberna@studenti.polito.it        %%%
%%%                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all
clear 
clc

tic;


%% Part 1 - Input data


[materials, sections, nodes, elements, restraints, forces] = txtimport();

%% Part 2 - Charateristic of each element

[elements] = elementsProperty(elements,nodes, sections, materials);

%% Part 3 - Original structure plot break

orig_structural_plotter(nodes, restraints, forces, elements)
uiwait


%% Part 4 - Structural solver

[nodes, d, restraints, forces, elements, sigma, disp, epsilon] = DisplmethSolver(materials, sections, nodes, elements, restraints, forces);

%% Part 5 - Plot

structural_plotter(nodes, d, restraints, forces, elements, sigma, disp, epsilon)

toc;




