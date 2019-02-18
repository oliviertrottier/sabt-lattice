function [treeout, occout, MaxT, timeseries] = sabt(alpha, beta, gamma, varargin)
% Function that simulates the growth of a dendritic tree on a triangular lattice.

% Branches grow at a rate of 1 lattice spacing/timestep.

% New branches are created uniformly across the tree and the number of new
% branches is determined by drawing from a Poisson distribution with rate =
% gamma (branching frequency).

% The persistence length of each branch is determined by beta. For each
% spatial step, the angle of growth changes by delta_theta where
% delta_theta is normally distributed with mean = 0 and variance =
% d/beta/2.
 
% Upon collision with another branch, branches retract a random length l, where
% l is exponentially distributed with mean = alpha. The retraction speed is equal to
% the growth speed. Retraction continues through the parent branch when 
% the sister branch has completely retracted.

% Stop on error for easy error handling
dbstop if error;

% Define the persistent variables.
%persistent d dy dx L N x y MaxL NBSoma MaxG MinB GrowthStartTime
persistent d dy dx L N x y
global tree occ
%% Options Handling
% Determine the input options with error handling when the input options
% are undefined.
options_name = {'Movie' 'MaxT', 'Timeseries'};
optionsvalue = zeros(1, numel(options_name));

input_options = varargin(cellfun(@isstr, varargin));
if ~ all(ismember(input_options, options_name))
    error('Input options are undefined.');
end

N_options = numel(options_name);
for i = 1:N_options
    [existtag, opt_loc] = ismember(options_name(i), cellfun(@num2str, varargin, 'Uni', 0));
    if existtag
        optionsvalue(i) = varargin{opt_loc + 1};
    end
end

% Construct the options structure.
options = struct();
for i = 1:N_options
    options.(options_name{i}) = optionsvalue(i);
end
%%
% Set the maximal number of time steps.
if options.MaxT > 0
    MaxT = options.MaxT;
else
    MaxT = 1200;
end

% Determine the number of branches at the soma.
NBSoma = 3;

% Determine the number of new growths at each time.
growthspertime = random('poisson', gamma, 1, MaxT - 1);

% Define the maximal number of branches and growths.
BranchingStartTime = 5;
MaxB = sum(2 * growthspertime(BranchingStartTime:end)) + NBSoma;
MaxG = ceil(MaxT * MaxB);
MaxL = round(500 / gamma);

% Construct the triangular lattice.
% The y value decreases with increasing row number.
if isempty(L) || L ~= MaxT + 100
    % Make the size of the lattice large enough so no branches hit the
    % boundary.
    L = MaxT + 100;
    d = 1;
    dy = d / 2;
    dx = sqrt(3) / 2 * d;
    N = round(2 * L / dy);
    [x, y] = meshgrid(- N / 2:1:N / 2, N / 2: - 1: - N / 2);
    x = x .* dx;
    y = y .* dy;
end

% Check if the lattice is smaller than the largest int16 integer.
if any(size(x) > 2 ^ 15)
    error('Error! The size of the lattice is bigger than 2^15.')
end

% Create a folder where the outputs will be saved.
savingfolder = cell2mat(regexp(mfilename('fullpath'), '/.*/', 'match'));
if isdir(savingfolder) == 0
    mkdir(savingfolder);
end

% Sample the retraction lengths and branching frequency and delay time.
retractionlengths = round(random('exp', alpha, 1, MaxB));
delta_theta = sqrt(d / (beta / 2)) * randn(MaxG, 1);

% Initialize the matrix that records the collision's position.
CollisionsPos = zeros(MaxB, 3);

% Define an array that will track the occupation of each lattice point.
% 0 corresponds to a free site.
% >=1 corresponds to an occupied site. The number indicates the branch ID
% that the site belongs to.
% BPocc defines the occupaction ID of a branch point.
BPocc = - 1;
occ = sparse(ceil((N + 1) / 2), ceil((N + 1) / 2), BPocc, N + 1, N + 1);

% Define the neighbours' search range in the x and y directions.
xupdrng = 2;
yupdrng = 2;

% Define the neighbours' vectors on a triangular lattice
neighbours_vec = int16([-1 1; -2 0; -1 -1; 1 -1; 2 0; 1 1]);
neighbours_ind = int16(sub2ind(2 * [xupdrng yupdrng] + 1, neighbours_vec(:, 1) + xupdrng + 1, neighbours_vec(:, 2) + yupdrng + 1));

% Find the boundary of the lattice.
rowmin = ceil((N + 1) / 2 - L / dx) + xupdrng;
rowmax = floor((N + 1) / 2 + L / dx) - xupdrng;
colmin = ceil((N + 1) / 2 - L / dy) + yupdrng;
colmax = floor((N + 1) / 2 + L / dy) - yupdrng;

% Initialize the number of growing branches.
NGrowths = NBSoma;
GrowingBranchesIDs = 1:NBSoma;
RetractingBranchesIDs = [];

% Define a structure to keep track of the tree branches.
tree = struct(...
    'PointsInd', cell(1, MaxB),... % Lattice indices of the points that form a branch.
    'TipPos', cell(1, MaxB),... % Position of the branch tip.
    'OriPos', cell(1, MaxB),... % Original Position of the branch tip.
    'ParentID', cell(1, MaxB),... % ID of the parent branch.
    'SisterID', cell(1, MaxB),... % ID of the sister branch.
    'DaughtersID', cell(1, MaxB),... % ID of the daugther branches.
    'State', num2cell(zeros(1, MaxB)),... % Growing (+1), shrinking (-1) or immobile (0).
    'Length', num2cell(zeros(1, MaxB, 'uint16')),... % Length of the branch.
    'Theta', cell(1, MaxB),... % Angle of growth.
    'Tau', cell(1, MaxB)); % Remaining amount of branch length to retract.

% Initialize the tree structure.
[tree.PointsInd] = deal(zeros(5 * MaxL, 2, 'int16'));
for i = 1:NBSoma
    tree(i).ParentID = 0;
    tree(i).SisterID = mod(i, 3) + 1;
    tree(i).OriPos = ceil((N + 1) / 2) .* ones(1, 2, 'int16');
    tree(i).TipPos = ceil((N + 1) / 2) .* ones(1, 2, 'int16');
    tree(i).PointsInd(1, :) = ceil((N + 1) / 2) .* ones(1, 2, 'int16');
    tree(i).State = 1;
end

% Define a structure to keep track of the timeseries data.
if options.Timeseries || options.Movie
    timeseries = struct();
    timeseries(MaxT) = struct();
else
    timeseries = [];
end

ll = 0;
mm = 0;
mm2 = 0;

%% Time loop
for t=1:MaxT
    % Randomize the order of the moving branches.
    MovingBranchesIDs = [GrowingBranchesIDs RetractingBranchesIDs];
    NMovingBranches = numel(MovingBranchesIDs);
    GrowthOrder = randperm(NMovingBranches)';

    for ii = 1:NMovingBranches
        % Get info about the current moving branch.
        BranchID = MovingBranchesIDs(GrowthOrder(ii));
        Length = tree(BranchID).Length;
        State = tree(BranchID).State;
        TipPos = tree(BranchID).TipPos;

        if State == 1
            % Find the occupation of neighboring lattice points, which may
            % influence the growth of the tip.
            minxInd = TipPos(2) - xupdrng;
            maxxInd = TipPos(2) + xupdrng;
            minyInd = TipPos(1) - yupdrng;
            maxyInd = TipPos(1) + yupdrng;
            EnvironmentOcc = occ(minyInd:maxyInd, minxInd:maxxInd);
            NeighOcc = EnvironmentOcc(neighbours_ind);

            % Check if the growing branch hits another branch
            if any(NeighOcc == 0)
                if Length < 1
                    % Select an empty site at random for the growth of the
                    % first point.
                    NeighSel = randsample(1:6, 1, true, full(NeighOcc == 0));
                    Theta = pi / 3 * (NeighSel - 1) + pi / 6;
                    occvalue = 0;
                else
                    % Determine the growth angle and round it to the
                    % closest integer multiple of pi/3 (60 deg).
                    ll = ll + 1;
                    Theta = tree(BranchID).Theta + delta_theta(ll);
                    NeighSel = round((mod(Theta, 2 * pi) - pi / 6) / (pi / 3) + 1);
                    occvalue = NeighOcc(NeighSel);
                end
            else
                % In this case, all neigbouring lattice points are occupied.
                occvalue = 1;
            end

            % Check for collision of the tip.
            if occvalue == 0
                % Find the new position of the tip.
                NewTipPos = TipPos + neighbours_vec(NeighSel, :);

                % Update the occupancy matrix.
                occ(NewTipPos(1), NewTipPos(2)) = BranchID;
                %occ(Neighbourslinind(NeighSel))=1;

                % Update the tree structure with the new tip position, branch length
                % and angle.
                tree(BranchID).TipPos = NewTipPos;
                tree(BranchID).Length = Length + 1;
                tree(BranchID).PointsInd(Length + 2, :) = NewTipPos;
                tree(BranchID).Theta = Theta;

            else
                % If occvalue=1, the growing tip has hit another branch.
                % Consequently, start the retraction of the tip.
                if alpha > 0
                    mm = mm + 1;
                    tree(BranchID).State = - 1;
                    tree(BranchID).Tau = retractionlengths(mm);

                    % Record the position of the collision.
                    CollisionsPos(mm - mm2, :) = [x(TipPos(1), TipPos(2)) y(TipPos(1), TipPos(2)) t];

                else
                    tree(BranchID).State = 0;
                end
            end

            % If BranchState~=1, then BranchState=-1. In this case, check
            % if the retraction is still ongoing (tau > 0).
        elseif tree(BranchID).Tau > 0
            % Update the occupancy and connectivity matrix.
            occ(TipPos(1), TipPos(2)) = 0;

            % Update the tree structure with the new tip position, branch
            % length and tau.
            tree(BranchID).Length = Length - 1;
            tree(BranchID).Tau = tree(BranchID).Tau - 1;
            tree(BranchID).TipPos = tree(BranchID).PointsInd(Length, :);
        else
            % If the retraction length reaches 0, stop the retraction.
            tree(BranchID).State = 0;
        end
    end

    % Next, handle the full retraction of moving branches.
    ZeroLengthIDs = MovingBranchesIDs([tree(MovingBranchesIDs).Length] == 0);
    rethandler(ZeroLengthIDs);

    % Store the timeseries data.
    if options.Timeseries
        % Find the lengths of the branches.
        Lengths = nonzeros([tree(1:NGrowths).Length]);
        NBranchPoints = nnz(occ == BPocc) - 1;

        % Calculate the radius of gyration.
        occupiedpos = occ ~= 0;
        Rg = rg([x(occupiedpos) y(occupiedpos)]);

        % Save other dynamics info.
        timeseries(t).Lengths = uint16(Lengths);
        timeseries(t).NGrowingTips = uint16(numel(GrowingBranchesIDs));
        timeseries(t).NRetractingTips = uint16(numel(RetractingBranchesIDs));
        timeseries(t).NCollisions = uint16(mm);
        timeseries(t).CollisionsPos = CollisionsPos(1:mm - mm2, :);
        timeseries(t).NBranchpoints = uint16(NBranchPoints);
        timeseries(t).Rg = Rg;

        mm2 = mm;
    end

    % Record Movie Data.
    if options.Movie ~= 0
        timeseries(t).tree = restructure(tree);
    end

    % When MaxT is reached, end the loop and do not perform branching for the
    % next timestep.
    if t >= MaxT
        break;
    end
    %% Branching.

    % When branching occurs, a random point on the tree is
    % selected as the branch point and a new moving tip is added at that point.
    TotalGrowths = growthspertime(t);

    if TotalGrowths > 0 && t >= BranchingStartTime
        [freeocc_rowInd, freeocc_colInd] = find(occ >= 1);
        NAvSites = size(freeocc_rowInd, 1);
        MaxGrowths = min(NAvSites, TotalGrowths);
        BranchingOrder = randsample(NAvSites, MaxGrowths);
        GrowthSitesInd = [freeocc_rowInd(BranchingOrder) freeocc_colInd(BranchingOrder)];

        for ii = 1:MaxGrowths
            % Find the position of the branch point, its
            % associated branch index and the position of the
            % branch point along this branch.
            BranchPointInd = int16(GrowthSitesInd(ii, 1:2));
            SeparatedBranchID = full(occ(BranchPointInd(1), BranchPointInd(2)));
            SeparatedBranchLength = tree(SeparatedBranchID).Length;
            SeparatedBranchPointsInd = tree(SeparatedBranchID).PointsInd(1:SeparatedBranchLength + 1, :);
            [~, BranchPointPos] = ismember(BranchPointInd,SeparatedBranchPointsInd,'rows');

            % If the branch point is located at the tip of one of the
            % branches, simply start the growth of that branch.
            if BranchPointPos == SeparatedBranchLength + 1
                tree(SeparatedBranchID).State = 1;
                continue;
            end

            % Add two new branches and growing tips.
            NGrowths = NGrowths + 2;
            NewBranch1ID = NGrowths - 1;
            NewBranch2ID = NGrowths;

            % Update the tree structure with the 1st new
            % branch. This branch corresponds to the downstream
            % segment of the separated branch.
            tree(NewBranch1ID).ParentID = SeparatedBranchID;
            tree(NewBranch1ID).SisterID = NewBranch2ID;
            tree(NewBranch1ID).DaughtersID = tree(SeparatedBranchID).DaughtersID;
            tree(NewBranch1ID).Theta = tree(SeparatedBranchID).Theta;
            tree(NewBranch1ID).Tau = tree(SeparatedBranchID).Tau;
            tree(NewBranch1ID).Length = SeparatedBranchLength - (BranchPointPos - 1);
            tree(NewBranch1ID).PointsInd(1:SeparatedBranchLength - (BranchPointPos - 1) + 1, :) = tree(SeparatedBranchID).PointsInd(BranchPointPos:SeparatedBranchLength + 1, :);
            tree(NewBranch1ID).State = tree(SeparatedBranchID).State;
            tree(NewBranch1ID).OriPos = BranchPointInd;
            tree(NewBranch1ID).TipPos = tree(SeparatedBranchID).TipPos;

            % Update the tree structure with the 2nd new
            % branch. This branch corresponds to the additional
            % segment.
            tree(NewBranch2ID).ParentID = SeparatedBranchID;
            tree(NewBranch2ID).SisterID = NewBranch1ID;
            tree(NewBranch2ID).Length = 0;
            tree(NewBranch2ID).PointsInd(1, :) = BranchPointInd;
            tree(NewBranch2ID).State = 1;
            tree(NewBranch2ID).OriPos = BranchPointInd;
            tree(NewBranch2ID).TipPos = BranchPointInd;

            % Update the separated branch's daughters IDs.
            DaughtersIndex = tree(SeparatedBranchID).DaughtersID;
            if ~ isempty(DaughtersIndex)
                for iii = 1:numel(DaughtersIndex)
                    tree(DaughtersIndex(iii)).ParentID = NewBranch1ID;
                end
            end

            % Update info of the separated branch.
            tree(SeparatedBranchID).Length = BranchPointPos - 1;
            tree(SeparatedBranchID).DaughtersID = [NewBranch1ID NewBranch2ID];
            tree(SeparatedBranchID).State = 0;
            tree(SeparatedBranchID).TipPos = BranchPointInd;

            % Update the occupancy matrix with the position of the branch
            % point.
            occ(BranchPointInd(1), BranchPointInd(2)) = BPocc;

            % Update the occupancy matrix where the core segments of
            % NewBranch1ID are located.
            NewBranch1IDLength = tree(NewBranch1ID).Length;
            for j = 2:NewBranch1IDLength
                PointInd = tree(NewBranch1ID).PointsInd(j, :);
                occ(PointInd(1), PointInd(2)) = NewBranch1ID;
            end

            % If the new branch has no daughters, its last occupancy must
            % also be updated.
            if isempty(tree(NewBranch1ID).DaughtersID) && NewBranch1IDLength > 0
                PointInd = tree(NewBranch1ID).PointsInd(tree(NewBranch1ID).Length + 1, :);
                occ(PointInd(1), PointInd(2)) = NewBranch1ID;
            end
        end

    end

    % Find the new set of branches that are growing and retracting.
    GrowingBranchesIDs = find([tree(1:NGrowths).State] == 1);
    RetractingBranchesIDs = find([tree(1:NGrowths).State] == - 1);
end

% Clean the tree structure by removing all non-necessary fields. Also,
% change the format of the recorde data to reduce memory needs.
treeout = restructure(tree);

% Change the occupancy matrix to a binary format (0=unoccupied, 1=occupied)
occout = logical(occ);

% Clear the global variables.
clear global tree occ;
end
