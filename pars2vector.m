function [param_vec, param_names] = pars2vector(pars, newpars)
% read in parameters values and 
% put into vector format
% Input: 
%    pars -- structure with parameters from set_params()
%    newpars -- binary, set to 1 if new parameters for printout
% Output:
%    parsvec -- vector of parameters values
%    param_names -- parameter names corresponding to vector values

param_names = fieldnames(pars);
numpars = length(param_names);
param_vec = zeros(size(param_names));
for ii = 1:numpars
    temp = strcat('pars.', param_names{ii});
    param_vec(ii) = eval(temp);
end

% optional: easier for writing up model equations and checking params
if newpars
    fid = fopen('newparamnames.txt','wt');
    for jj = 1:numpars
        fprintf(fid, '%s = params(%d);\n',param_names{jj},jj);
    end
    fprintf('TO DO: copy results from newparamnames.txt into parameters set up in model function \n')
end
fclose('all');
end

