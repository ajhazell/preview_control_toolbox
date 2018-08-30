function Value = get(sys,Property)
%standard 'get' data accessor


switch Property
	case 'n'
		Value=sys.n;
	case 'p'
		Value=sys.p;
	case 'q'
		Value=sys.q;
	case 'l'
		Value=sys.l;
	case 'm'
		Value=sys.m;
    case 'B1'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);% This is lazy - should really compute B1 from scratch
		Value=B1;
    case 'B2'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=B2;
    case 'C1'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=C1;
    case 'C2'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=C2;
    case 'D11'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=D11;
    case 'D12'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=D12;
    case 'D21'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=D21;
    case 'D22'
        [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
		Value=D22;

    otherwise
		Value=get(sys.ss,Property); % pass it up the class structure
end
