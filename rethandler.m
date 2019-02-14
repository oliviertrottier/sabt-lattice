function rethandler11(MovingBranchesIDs)
% Function that handles the full retraction of branches.
global tree
global occ

NBranches = numel(MovingBranchesIDs);
for i = 1:NBranches
    % Get info on the retracting branch.
    BranchID = MovingBranchesIDs(i);
    SisID = tree(BranchID).SisterID;
    ParID = tree(BranchID).ParentID;
    SisBranchLength = tree(SisID).Length;
    
    % If the branch is at the soma, simply stop its growth. Soma branches
    % have a parent ID of 0.
    if ParID == 0
        tree(BranchID).State = 0;
        continue;
    end
    
    if SisBranchLength > 0
        % Add the length and branch points of the sister to the parent branch.
        ParBranchLength = tree(ParID).Length;
        tree(ParID).Length = ParBranchLength + SisBranchLength;
        tree(ParID).PointsInd(ParBranchLength+1:tree(ParID).Length+1, :) = tree(SisID).PointsInd(1:SisBranchLength+1, :);
        tree(ParID).TipPos = tree(SisID).TipPos;
        tree(ParID).DaughtersID = [tree(SisID).DaughtersID];
        tree(ParID).Tau = tree(SisID).Tau;
        tree(ParID).Theta = tree(SisID).Theta;
        tree(ParID).GrowthDelay = tree(SisID).GrowthDelay;
        tree(ParID).State = tree(SisID).State;
        
        % Update the parent ID of the daughters.
        DaughtersIndex = [tree(SisID).DaughtersID];
        NDaughters = numel(DaughtersIndex);
        for j = 1:NDaughters
            tree(DaughtersIndex(j)).ParentID = ParID;
        end
        
        % Update the occupancy matrix where the core segments of the sister
        % branch are located.
        for j = 1:SisBranchLength
            PointInd = tree(SisID).PointsInd(j, :);
            occ(PointInd(1), PointInd(2)) = ParID;
        end
        
        % If the sister branch has no daughters, its last occupancy must
        % also be changed.
        if NDaughters == 0
            PointInd = tree(SisID).PointsInd(SisBranchLength+1, :);
            occ(PointInd(1), PointInd(2)) = ParID;
        end
        
        % Set the length of the sister branch to zero.
        tree(SisID).Length = 0;
    elseif BranchID < SisID
        % In this case, both branches retracted at the same time. In
        % such case, the retraction of the parent branch starts and its
        % retraction length is given by the largest retraction length
        % of its daughters.
        %         if tree(ParID).Length==0
        %             error('error')
        %         end
        ParTau = max(tree(BranchID).Tau, tree(SisID).Tau);
        if ParTau > 0
            tree(ParID).Tau = ParTau;
            tree(ParID).State = -1;
        end
        
        % Remove the daughters of the parent.
        tree(ParID).DaughtersID = [];
        
        % Change the occupancy of the original branch point back to a
        % normal occupied site (1).
        PointInd = tree(BranchID).OriPos;
        occ(PointInd(1), PointInd(2)) = ParID;
    end
    
    % Stop the motion of the branch and its sister.
    tree(BranchID).State = 0;
    tree(SisID).State = 0;
    
end
end