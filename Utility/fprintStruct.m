function [] = fprintStruct(fileID, mystruct)
%
%fprintStruct prints a structure (or a vector of structures) to a text file.
%   The structure can contain fields of different types, including vectors.
%   Use fopen(FILENAME) to get fileID.
%
%   [] = fprintStruct(fileID, mystruct)
%
%   Inputs:
%       fileID      ID of the destination text file (integer)
%       mystruct    Struct to print to the file     (struct)
%
%   See also: fopen, fprintf
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

arguments
    fileID      double  {mustBeInteger} 
    mystruct    struct
end

if length(mystruct) > 1 && ~isnumeric(mystruct)
    for i=1:length(mystruct)
        fprintStruct(fileID,mystruct(i))
    end
else
    
    fn = fieldnames(mystruct);
    for k=1:numel(fn)
        if( length(mystruct.(fn{k})) > 1 )
            fprintf(fileID,'%s= %s\n',fn{k},mat2str(reshape(mystruct.(fn{k}), 1, []) ,3));
            
        elseif( isinteger(mystruct.(fn{k})) )
            fprintf(fileID,'%s= %d\n',fn{k},mystruct.(fn{k}));
            
        elseif( isnumeric(mystruct.(fn{k})) )
            fprintf(fileID,'%s= %.2f\n',fn{k},mystruct.(fn{k}));
            
        elseif( islogical(mystruct.(fn{k})) ) 
            fprintf(fileID,'%s= %s\n',fn{k},mat2str(mystruct.(fn{k})));
            
        elseif( isa(mystruct.(fn{k}),'function_handle') )
            fprintf(fileID,'%s= %s\n',fn{k},func2str(mystruct.(fn{k})));
            
        elseif( iscell(mystruct.(fn{k})) )
            fprintf(fileID, '%s= %s', fn{k},mat2str(mystruct.(fn{k}){:},3));
            fprintf(fileID, '\n');
            
        else
            fprintf(fileID,'%s= %s\n',fn{k},mystruct.(fn{k}));
        end
        
    end
    fprintf(fileID,'\n');
    
end

end

