% Example script to run the sabt script that creates a dendritic arbor
%% Run the simulation
rng(1);
alpha = 50;
beta = 100;
gamma = 1;
[Tree,occ,TotalTimeFrames,timeseries] = sabt(alpha,beta,gamma);

%% Plot the tree
% Form the triangular lattice
d=0.4; % lattice spacing (0.4 um)
dy=d/2;
dx=sqrt(3)/2*d;
N=size(occ,1)-1;
[x,y]=meshgrid(-N/2:1:N/2,N/2:-1:-N/2);
x=x.*dx;
y=y.*dy;

% Initialize the plot.
figure
a=gca;
xlabel('x (\mum)')
ylabel('y (\mum)')

% Define a line plot for each branch.
NBranches = numel(Tree);
for i=1:NBranches
    BranchLength = Tree(i).Length;
    LatticeInd = Tree(i).PointsInd(1:BranchLength+1,:);
    line('Parent',a,'XData',x(1,LatticeInd(:,2))','YData',y(LatticeInd(:,1),1))
end
max_x = max(abs([a.Children.XData]));
max_y = max(abs([a.Children.YData]));
max_L = max(max_x,max_y);
xlim(max_L*[-1 1])
ylim(max_L*[-1 1])