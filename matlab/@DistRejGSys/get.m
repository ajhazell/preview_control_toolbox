function Value = get(sys,Property)
% standard 'get' accessor function

switch Property
	case 'lr'
		Value=sys.lr;
	case 'lw'
		Value=sys.lw;
	otherwise
		Value=get(sys.GenSys,Property); % pass it up the class structure
end
