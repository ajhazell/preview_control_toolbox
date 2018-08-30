function result = subsref(sys,Struct)
%SUBSREF Subscripted reference for all props with 'get' defined
%

ni = nargin;
if ni==1,
   result = sys;
   return
end
StructL = length(Struct);

% Peel off first layer of subreferencing
switch Struct(1).type
case '.'
   % The first subreference is of the form sys.fieldname
   % The output is a piece of one of the system properties
   try
      if StructL==1,
         result = get(sys,Struct(1).subs);
      else
         result = subsref(get(sys,Struct(1).subs),Struct(2:end));
      end
   catch
      rethrow(lasterror)
   end
otherwise
   error(sprintf('Unknown reference type: %s',Struct(1).type));
end
