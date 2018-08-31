# preview_control_toolbox
Matlab Preview Control Toolbox v0.3
===================================

This is the first new version for 10 years! The core algorithms are unchanged, but it is now updated to work with the modern Matlab class syntax and tested against Matlab 2018a. 

Full technical details may be found in: 

An efficient algorithm for discrete-time Hinfinity preview control, A. Hazell and D. Limebeer, Submitted to Automatica.

and:

A design framework for discrete-time H2 preview control, A. Hazell and D. Limebeer, Submitted to Trans. ASME, Journal of Dynamic Systems, Measurement and Control.

Notation and signal names follow that used in these documents.

Classes included in this version are

GenSys		: Generalised plant
DistRejGSys     : Generalised plant with exogenous input split into w and r components
DistRejPrevSys  : Generalised plant for the general previewable disturbance rejection problem
PrevTrackSys    : Generalised plant for the preview tracking problem

The class hierachy is

		GenSys
		   |		
	      DistRejGSys
		   |    
	    DistRejPrevSys
		   | 
	      PrevTrackSys
		   |
              LQRTrackSys   


Please add an 'issue' to the issue tracker if you have any comments or suggestions.


Getting Started
===============

Add the root folder 'matlab' to your MATLAB path, and see 'PCTDemo.m' for an introduction to using this code.

Online help can can be obtained by, for example, typing 
	>>help PrevTrackSys

Further documentation will follow in the next version. However, I've been saying that for the last 10 years and haven't done it. If you've got any specific questions please add an issue to the issues list and I'll try to get back to you.


