function matrix=rep(matrix,oldvalues,newvalues)
if ~isempty(matrix)
    [~,locs]=ismember(matrix,oldvalues);
    nonzerolocs=locs~=0;
    matrix(nonzerolocs)=newvalues(locs(nonzerolocs));
    %
    %     replacements=find(ismember(oldvalues,matrix));
    %     N_replacements=numel(replacements);
    %     if N_replacements>0
    %         matrix_temp=matrix;
    %         for i=1:N_replacements
    %             matrix(matrix_temp==oldvalues(replacements(i)))=newvalues(replacements(i));
    %         end
    %     end
end
end