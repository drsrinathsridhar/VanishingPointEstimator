Description
===============
VanPoints is a class that implements the computation of vanishing points from a
given image using the approach proposed in the research paper 'Non-Iterative
Approach for Fast and Accurate Vanishing m_Point Detection' by Jean-Philippe
Tardif. The paramers for the algorithm are the minimum line length and J-Linkage
model size. Please see the original paper for explanantion on these. The methods
in the VanPoints class are self-explanatory. The primary method is the
findVanishingPoints(). Please refer to sample_usage.cpp for help on using the
class.


Building
===============
Requires OpenCV (2.1 and 2.2 tested). On linux systems please do the following

``make``
``./findVanPoints threepp.jpg 10 500``

Also tested on VS2008 on a Windows XP 64 bit machine. Please don't forget to add
the OpenCV include paths and the OpenCV libraries to your project properties.
Please email Srinath in case of problems.


Revision History
==============================

Version 1.0: Original written by Yu Xiang and Srinath Sridhar for the course
						 EECS 442, Fall 2010 at the University of Michigan, Ann Arbor
Version 2.0: For private use
Version 3.0: Converted to a class and open sourced by Srinath

Contact: srinaths@umich.edu
WWW: umich.edu/~srinaths

