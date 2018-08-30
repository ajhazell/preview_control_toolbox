# preview_control_toolbox
Matlab Preview Control Toolbox 
==========================

This is the first release of a toolbox for Preview Control.

This toolbox is not yet feature complete, and whilst it has been thoroughly tested for accuracy of the algorithms, it has not yet undergone testing for software robustness.

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
PrevTrackSys    : Generalised plant for the LQR preview tracking problem

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


Please direct any comments/suggestions to a.hazell@ic.ac.uk


Getting Started
===============

Add the root folder 'PCT' to your MATLAB path, and see 'PCTDemo.m' for an introduction to using this code.

Online help can can be obtained by, for example, typing 
	>>help PrevTrackSys

Further documentation will follow in the next version.


