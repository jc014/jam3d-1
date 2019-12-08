Simulator
=========

Examples of how to perform a simulation is provided in the folder ::

  jam3d/simulation

workflow
--------

1) Use the template.xlsx to define the kinematics and setup the relative uncertainties 
2) Run the script simulation.py (./simulation.py). It will output a file simulation.xlsx

comments
--------

- The script simulation.py is generic and can be modified easily. It should be viewed as 
  an example. So the naming for the input  and output xlsx files can be changed if desired. 
  This is useful if there are more than one kind of data sets to be simulated

- In the example the simulated data coincides with the theory. One can add gaussian noise 
  if desired via ::

    alpha2=simdata['stat_u']**2 + simdata['syst_u']**2
    alpha=alpha2**0.5
    simdata['value']+=alpha*np.random(len(alpha)) 
  
available examples
------------------

- example00: unpolarized FUU in SIDIS



