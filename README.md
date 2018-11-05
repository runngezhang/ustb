# UltraSound ToolBox (USTB) #

An open source MATLAB toolbox for beamforming, processing, and visualization of ultrasonic signals. The USTB is developed as a joint effort of:
 
* [Department of Circulation and Medical Imaging of NTNU](https://www.ntnu.no/isb), 
* [Department of Informatics of the University of Oslo](http://www.uio.no/), and
* [CREATIS Laboratory of the University of Lyon](https://www.creatis.insa-lyon.fr/site7/en).

### How do I get set up? ###

* Just clone the repository and add the folder (without subfolders) to MATLAB's path

### Current version ###

The USTB is still under development, so there might be larger structural changes. The current version in main is;

* v2.1: https://bitbucket.org/ustb/ustb/commits/tag/v2.1.2

compared to the previous version

* v2.0: https://bitbucket.org/ustb/ustb/commits/tag/v2.0.1

The main changes are:

* updated the processing pipeline defining the pre, mid and postprocess objects 
* the apodization object have been rewritten
* a GPU implementation of the DAS beamformer have been added 
* several bugfixes and other improvements have been done as well.

### Documentation ###
Unfortunately, we have not had the time or resources to write a full documentation of the USTB. However, there are plenty of well documented examples that will help you to get started and hopefully understand the code.

### How to contribute? ###
First off all, please make your self familiar with the Gitflow workflow. See for example this tutorial: https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow. According to the Gitflow workflow there are two types of contributions to the code a __hotifx__ and a __feature__.

* __hotfix__: This is a change done on the main branch. Typically, an urgent bugfix often related to a reported issue: https://bitbucket.org/ustb/ustb/issues?status=new&status=open.
* __feature__: This s a change done on the develop branch. Typically, a larger change to the code or added or improved functionality. 

To contribute to the project with your code you should do the following:

* Step 1: Create your own fork of the project. 
* Step 2: Create a hotfix branch from the main branch to fix an urgent bug, or a feature branch from the develop branch to add a feature.
* Step 3: Create a pull request from your forked repository back to our repository and add @omrindal and @alfonsomolares as reviewers. 
* Step 4: Once we have time to review your code and are happy with the changes we will merge your pull request into the USTB and you have sucessfully contributed!

Once we are happy and confortable that the develop branch is stable and useful we will merge it into main and a new version will be released :D

### Who do I talk to? ###

The project administrators are:

* Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>,
* Ole Marius Hoel Rindal <omrindal@ifi.uio.no>
* Olivier Bernard <olivier.bernard@insa-lyon.fr> 
 

Collaborators:

* Andreas Austeng 
* Arun Nair
* Muyinatu A. Lediju Bell, 
* Lasse Løvstakken 
* Svein Bøe 
* Hervé Liebgott 
* Øyvind Krøvel-Velle Standal 
* Jochen Rau 
* Stefano Fiorentini
