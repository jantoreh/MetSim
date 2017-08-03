% The MetSim program
% Author: Jan-Tore Horn
% Date: August 2017
% License: GPL

function X = metsim(varargin)


% Sampling from an environmental distribution given the dependency table in
% Horn (2018) (to be submitted).
% 

% Read inputfile
inpfile = 'metsim.inp';
fid = fopen(inpfile,'r');

fgetl(fid); % Dummy
n = fscanf(fid,'%f',1);fgets(fid); % Number of variables
varnames = textscan(fid,repmat('%s ',1,n),1,'CommentStyle','-'); % Variable names
varunits = textscan(fid,repmat('%s ',1,n),1,'CommentStyle','-'); % Variable units
vardists = textscan(fid,repmat('%s ',1,n),1,'CommentStyle','-'); % Variable distributions
depfile = fscanf(fid,'%s ',1);fgetl(fid); % Dependency file
paramfile = fscanf(fid,'%s ',1);fgetl(fid); % Parameter file
fgetl(fid); % Dummy
N_init = fscanf(fid,'%d',1);fgetl(fid); % Rejected samples
init = fscanf(fid,'%f',n);fgetl(fid); % Starting point
N = fscanf(fid,'%d',1);fgetl(fid); % Output results
lower = fscanf(fid,'%f',n);fgetl(fid); % Lower limit
upper = fscanf(fid,'%f',n);fgetl(fid); % Upper limit
seed = fscanf(fid,'%d',1);fgetl(fid); % Seed number
fgetl(fid);
outputname = fscanf(fid,'%s',1);fgetl(fid);
plotdepfile = fscanf(fid,'%s',1);fgetl(fid);
fclose(fid);

% Default values
plott = 1;

% Manual inputs overriding inputfiles
switch nargin
    case 1
        N = varargin{1};
    case 2
        N = varargin{1};
        plott = 0;
end

% Make Octave find functions
isoctave = (exist('OCTAVE_VERSION')~=0); % Checks if octave running
if ~isoctave
    rng(seed);
else
    rand ('seed', 200);
    addpath('./private')
    pkg load statistics % Does not contain gmdistribution yet..
    error('Not fully compatible with Octave yet..')
end


% Begin
dep = load(depfile);
param = load(paramfile);
plotdep=load(plotdepfile);
for i =1:n  % Loop over all variables
    numpar(i) = type2numpar(vardists{i}{1});
end

% Check for errors
if any(size(dep)~=n)
    error('Invalid number of dependencies in dependency file')
end
if any(sum(dep,2)>2)
    error('Too many dependencies for some variables. Maximum is 1.')
end
if size(param,2)~=7
    error('The parameter file must have 7 columns')
elseif size(param,1)~=sum(numpar)
    error(['The current setup requires ',num2str(sum(numpar)),' rows in the parameter file'])
elseif any(triu(dep)-eye(size(dep,1))>0)
    error('Upper triangle in dependency matrix must be zero')
end

% Remove autocorrelation
dep = dep - eye(size(dep));

if any(init<lower) % Avoid starting outside wanted interval
    init(init<lower)=lower(init<lower)+0.001;
end

X = zeros(N+N_init,n); % Pre-allocate
x = init; % Starting values
i = 1;
nf = 1000; % Number of points to evaluate PDF

while i<(N+N_init+1)
    X(i,:) = x;
    i = i+1;
    
    % Sample
    e = rand(n,1);
    
    % Loop over variables
    for j=1:n
    
            [~,depvar]=max(dep(j,:)); % Index of dependency
            p = zeros(numpar(j),1);
            id1 = sum(numpar(1:j))-numpar(j)+1; % Row index start
            id2 = sum(numpar(1:j)); % Row index stop

         if sum(dep(j,:)) > 0 % Dependent variable
            p=paramfit(x(depvar),param(id1:id2,:));
         else
             p=param(id1:id2,1);
         end
         
         % Get PDF
         vals = linspace(lower(j),upper(j),nf);
         f = density(vardists{j}{1},vals,p);

         % Create CDF with unique insertions
         F = cumsum(f);
         F(1)=0;
         F=F/F(end);
         [F,index] = unique(F);
         vals = vals(index);
         
            
         % Sample new x
         x(j) = interp1(F,vals,e(j));
         
    end
    
    
    
end
X(1:N_init,:)=[]; % Remove initialization



%% Post-processing

% Save seastates
dlmwrite([outputname,'.out'],X,'delimiter',' ','precision','%10.2f');

% Exit if plotting is not requested
if plott ~= 1
    return
end

%% Create animation
fig=figure;
set(fig,'position',[100,20,800,600]);
gifname=[outputname,'.gif'];
color = [0,0,1];
grid on
dx = 20;
first=1;

% number of subplots
nx = 3;
nplot = sum(sum(plotdep));
ny = ceil(nplot/nx);

[idx,idy]=find(plotdep>0);

for i=2:dx:size(X,1)
   
    for j=1:nplot
    
    subplot(nx,ny,j)
    scatter(X(1:i,idx(j)),X(1:i,idy(j)),'.','MarkerEdgeColor',color);
    axis([lower(idx(j)),upper(idx(j)),lower(idy(j)),upper(idy(j))]);hold on
    scatter(X(i,idx(j)),X(i,idy(j)),'*','MarkerEdgeColor',[1,0,0])
    x = linspace(lower(idx(j)),upper(idx(j)),100);
    dx = max(x)-min(x);
    y = linspace(lower(idy(j)),upper(idy(j)),100);
    dy = max(y)-min(y);
    px=ksdensity(X(1:i,idx(j)),x);
    py=ksdensity(X(1:i,idy(j)),y);
    plot(x,px/max(px)*0.2*dy+min(y),'color','green','linewidth',1)
    plot(py/max(py)*0.2*dx+min(x),y,'color','green','linewidth',1)
    ylabel([varnames{idy(j)}{1},' [',varunits{idy(j)}{1},']'])
    xlabel([varnames{idx(j)}{1},' [',varunits{idx(j)}{1},']'])
    hold off
    grid on;box on
    
    end
    
    drawnow
    
    % Capture
    frame=getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256); 
    
     % Write to the GIF File 
      if first==1
          delete(gifname)
          imwrite(imind,cm,gifname,'gif', 'Loopcount',1); 
          first=0;
      else 
          imwrite(imind,cm,gifname,'gif','WriteMode','append'); 
      end 
    
end








