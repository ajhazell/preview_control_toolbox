function G1G2=mtimes(G1,G2)
% G3=mtimes(G1,G2)
% 
% Overloading of the multiplication operator for 
% GenSys*ss or ss*GenSys


if isa(G1,'GenSys') && isa(G2,'ss')
	if isa(G2,'GenSys')
		error('Multiplication of GenSys*GenSys is not defined')
	end
	[n1,p1,q1,l1,m1]=Getsz(G1);
	G1G2=GenSys(G1.ss*G2,q1,m1);
	
elseif isa(G2,'GenSys') && isa(G1,'ss')
	if isa(G1,'GenSys')
		error('Multiplication of GenSys*GenSys is not defined')
	end
	[n2,p2,q2,l2,m2]=Getsz(G2);
	G1G2=GenSys(G1*G2.ss,q2,m2);
else
	error('Incompatible classes')
end
	
	
