Tutorial
========

dependencies
------------

- Only Linux and OSX is supported

- We recommend to install anaconda (python2) which 
  comes will all the necessary libraries

installation
------------

Clone the codes from https://github.com/JeffersonLab/jam3d

Certain enviroment variables needs to be set in your terminal session. 
For bash they are ::

  export JAM3D=path2fitpack
  PYTHONPATH=$JAM3D:$PYTHONPATH
  export PATH=$JAM3D/bin:$PATH

Alternatively you can source the setup files ::

  source  setup.bash 

Once the enviroment variables are set the scripts can be called from 
anywhere in your system. In other words, there is no need to work within 
the same code folder. As a good practice, create dedicated folders for a
given analysis. 

how to use the code
-------------------

There are two ways to use the code:

1) As a TMD extraction tool

2) As a simulator for TMD observables



