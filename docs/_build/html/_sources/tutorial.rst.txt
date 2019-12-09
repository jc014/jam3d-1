Tutorial
========

dependencies
------------

- Only Linux and OSX is supported

- We recommend to install anaconda (python2) which 
  comes will all the necessary libraries

clonning the code
-----------------

Go to https://github.com/JeffersonLab/fitpack
and fork the repository under your GitHub account 
and clone your forked fitpack repo 

To sync your fork with updates on upstream we need to first 
setup the upstream ::

  git remote add upstream git@github.com:JeffersonLab/jam3d.git
  git fetch upstream

After that you can use ::
  
  git pull upstream master

to be able to sync your fork.


installation
------------

Go to  https://github.com/JeffersonLab/jam3d
and fork the repository under your GitHub account 
and clone your forked fitpack repo 

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


How to contribute
-----------------

You should run, test and develop the code under your fork. 
If you want to share some develpment, you can send us a pull 
request






