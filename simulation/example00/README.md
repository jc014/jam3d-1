# simulation for unpolarized SIDIS

## getting started

- Use the template.xlsx to define the kinematics and setup the relative uncertainties 
- Run the script simulation.py. It will output a file simulation.xlsx
- Add the dataset into input.py 

## additional comments

- The script simulation.py is generic and can be modified easily. It should be viewed as 
  an example. So the naming for the input  and output xlsx files can be changed if desired. 
  This is useful if there are more than one kind of data sets to be simulated

- In the example the simulated data coincides with the theory. One can add gaussian noise 
  if desired via:

  ```
  alpha2=simdata['stat_u']**2 + simdata['syst_u']**2
  alpha=alpha2**0.5
  simdata['value']+=alpha*np.random(len(alpha)) 
  ```
  
  

