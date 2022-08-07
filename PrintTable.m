function fid = PrintTable(sys)
    fid = fopen('table.txt','w');
    %
    % print the variable names
    %
    fprintf('           ');
    nc = length(sys.variables);
    for i = 1:nc;
        fprintf('%5s ',sys.variables{i});
    end
    fprintf('\n');
    %
    % print the variable numbers
    %
    fprintf('           ');
    for i = 1:nc
        fprintf('%5i ',i);   
    end
    fprintf('\n');
    % print a row separator
    fprintf('           ');
    for i = 1:nc
        fprintf('------');
    end
    fprintf('\n');
    
    % print the total mass stoichiometric coefficients
    M = sys.M;
    nr = size(M,1);
    for j = 1:nr
        fprintf('%11s',sys.mass{j});
        for i = 1:nc
            fprintf('%5i ',M(j,i));
        end
        fprintf('\n');
    end
    fprintf('           ');
    % print a row separator
    for i = 1:nc
        fprintf('------');
    end
    fprintf('\n');
    % print the acid-base coefficients
    K = sys.K;
    %nr = size(K,1);
    nr = size(sys.system);
    for j = 1:nr
        fprintf('%11s',sys.system{j});
        for i = 1:nc
            fprintf('%5i ',K(j,i));
        end
        fprintf('\n');
    end
    fprintf('\n');
    fclose(fid);
end
