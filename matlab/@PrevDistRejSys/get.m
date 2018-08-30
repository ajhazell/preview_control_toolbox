function Value = get(sys,Property)
% provides standard 'get' data access function

%pnames not used because I want these values to be read only

switch Property
	case 'N'
		Value=sys.N;
	case 'lw'
		Value=sys.lw;
	case 'lr'
		Value=sys.lr;
	case 'nwr'
		Value=sys.nwr;
	case 'ng'
		Value=sys.ng;
	case 'Wr'
		Value=sys.Wr;
	case 'Wz'
		Value=sys.Wz;
	case 'Ww'
		Value=sys.Ww;
	case 'G'
		Value=sys.G;
    case 'GW'
		Value=sys.GW;
    case 'qg'
        Value=sys.GW.q;
	otherwise
		Value=get(sys.DistRejGSys,Property); % pass it up the class structure
end
