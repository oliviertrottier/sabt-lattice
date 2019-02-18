function tree = restructure(treein)
% Function that restructures the input structure treein to remove
% zero-length branches and unnecessary fields.
tree = treein;
NBranches = numel(tree);

% Define the current IDs.
IDs = 1:NBranches;

% Find the number of branches at the soma and initialize their depth to 1.
BranchesID = find([tree.ParentID] == 0);
depth = 1;

% Arrange the depths and color of the branches.
while ~isempty(BranchesID)
    % Assign the depth.
    [tree(BranchesID).Depth] = deal(depth);
    
    % Find the branches of the next depth. These will be the daughters of
    % all current branches.
    BranchesID = [tree(BranchesID).DaughtersID];
    
    % Increment the depth.
    depth = depth + 1;
end

% Remove the branches of zero length and recount the number of branches.
Zerobranches_Ind = [tree.Length] == 0;
tree(Zerobranches_Ind) = [];
IDs(Zerobranches_Ind) = [];
NBranches = numel(tree);

% End the sorting if the tree is empty.
if NBranches == 0
    return;
end

% Reassign the sister IDs when the number of branches is less than 3.
if NBranches < 3
    switch NBranches
        case 2
            tree(1).SisterID = 2;
            tree(2).SisterID = 1;
        case 1
            tree(1).SisterID = 0;
    end
end

% Reorder the tree structure by ascending depth.
[~, NewOrder] = sort([tree.Depth]);
tree = tree(NewOrder);

% Replace the oldIDs by the newIDs.
oldIDs = IDs(NewOrder);
newIDs = 1:NBranches;
ParIDs = rep([tree.ParentID], oldIDs, newIDs);
SisIDs = rep([tree.SisterID], oldIDs, newIDs);
DaughtersIDs = rep([tree.DaughtersID], oldIDs, newIDs);
j = 1;

for i = 1:NBranches
    tree(i).ParentID = ParIDs(i);
    tree(i).SisterID = SisIDs(i);
    
    if ~isempty(tree(i).DaughtersID)
        tree(i).DaughtersID = DaughtersIDs(2*j-1:2*j);
        j = j + 1;
    end
end

% Check if all IDs have properly been replaced.
allIDs = cell2mat({[tree.ParentID], [tree.SisterID], [tree.DaughtersID]});
if any(allIDs > NBranches)
    disp(['# Branches = ', NBranches])
    disp(allIDs(allIDs > NBranches))
    error('Some ID replacements are incorrect.')
end

% Reduce memory usage.
MinMemFunc = @(x) str2func(['uint', num2str(2^(max(3, ceil(log2(log2(double(max(x))))))))]);

ID_mem_reduce = MinMemFunc(NBranches);
Depth_mem_reduce = MinMemFunc([tree.Depth]);
Length_mem_reduce = MinMemFunc([tree.Length]);

for i = 1:NBranches
    % Remove the non-necessary points indice
    tree(i).PointsInd = uint16(tree(i).PointsInd(1:tree(i).Length+1, :));
    
    % Dhange data type to reduce memory.
    tree(i).ParentID = ID_mem_reduce(tree(i).ParentID);
    tree(i).SisterID = ID_mem_reduce(tree(i).SisterID);
    tree(i).DaughtersID = ID_mem_reduce(tree(i).DaughtersID);
    tree(i).Depth = Depth_mem_reduce(tree(i).Depth);
    tree(i).Length = Length_mem_reduce(tree(i).Length);
end

% Remove unnecessary fields.
fields_tokeep = {'PointsInd', 'ParentID', 'SisterID', 'DaughtersID', 'Length', 'Depth', 'Color'};
tree_fieldnames = fieldnames(tree);
tree = rmfield(tree, tree_fieldnames(~ismember(tree_fieldnames, fields_tokeep)));
end