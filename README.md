    #########################
    ###                   ###
    ###    TrussSolver    ###
    ###                   ###
    #########################
    
:Author: AntonioSberna
:Date: 12/04/2019
:Revision: 0.1
:License: Public Domain

= Project:

It is a FEM solver that calculates displacements, stresses and reactions of a plane truss system made of isotropic material subjected to nodal forces.

== Step 1: Installation
Please describe the steps to install this project.

For example:

1. Open this file
2. Edit as you like
3. Release to the World!



=== Folder structure

....
 
  ├── txtimport.m                 =>  import data from txt files
  ├── elementsProperty.m          =>  elaborate the raw input of each bar
  ├── DisplmethSolver.m           =>  fem solver using direct stiffness matrix method
  ├── orig_structural_plotter.m   =>  it draws the non-deformed shape of the structures
  ├── structural_plotter.m        =>  it draws the deformed shape of the structures with stress and displacements

  └── txtexport.m                 =>  export all the outputs in txt files
....

=== License
This project is released under a {License} License.

=== Contributing
To contribute to this project please contact AntonioSberna 

 
To define the system you simply have to modify the six .txt files as shown in the first two rows of each document with material and sections property, nodes coordinates, connectivity table, restraints and at the end nodal forces applied (don't modify file names).

Then you have to run TrussSolver.m matlab script to check the setting and collect the results (as .txt output and graphically).
