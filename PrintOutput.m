function out = PrintOutput(varargin)
    if (nargin == 2);
        sys = varargin{1};
        fid = varargin{2};
        nv = length(sys.variables);
        fprintf(fid,'%15s ','              ');
        for i = 1:nv
            fprintf(fid,'%16s ',sys.variables{i});
        end
        fprintf(fid,'\n');
        fprintf(fid,'%15s ','  ');
        for i = 1:nv
            fprintf(fid,'%16s ','---------------');
        end
        fprintf(fid,'\n');
    else
        y    = varargin{1};
        sigy = varargin{2};
        yobs = varargin{3};
        wobs = varargin{4};
        sys  = varargin{5};
        fid  = varargin{6};
        p = sys.p;
        q = sys.q;
        nv = length(sys.variables);
        fprintf(fid,'%15s ','prior');
        for i = 1:nv
            fprintf(fid,'%7.4f ± %6.4f ',yobs(i),wobs(i).^(-0.5));
        end
        fprintf(fid,'\n');
        fprintf(fid,'%15s ','posterior');
        for i = 1:nv
            fprintf(fid,'%7.4f ± %6.4f ',y(i),sigy(i));
        end
        fprintf(fid,'\n');
        fprintf(fid,'%15s ','  ');
        for i = 1:nv
            fprintf(fid,'%16s ','---------------');
        end
        fprintf(fid,'\n');
    end
    out = [];
end

