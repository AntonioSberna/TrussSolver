
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


%% Part 1 - Input data

[materials, sections, nodes, elements, restraints, forces] = txtimport();

%% Part 2 - Structural solver

[nodes, d, restraints, forces, elements, sigma, disp, epsilon] = ...
    DisplmethSolver(materials, sections, nodes, elements, restraints, forces);

%% Part 3 - Plot

structural_plotter(nodes, d, restraints, forces, elements, sigma, disp, epsilon)





