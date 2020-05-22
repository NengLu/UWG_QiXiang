%-------------------------------------------------------------------------%
% References: 
% https://topotoolbox.wordpress.com/2015/01/22/plotting-colored-stream-networks/
% https://topotoolbox.wordpress.com/2017/09/21/an-alternative-to-along-river-swath-profiles/
%-------------------------------------------------------------------------%

DEM = GRIDobj('Data_DEM/srtm_TRR_90m_input.tif');
th_discharge = 50000; 

DEMf = fillsinks(DEM);
FD = FLOWobj(DEMf);

A = flowacc(FD);
W = A>th_discharge;
S = STREAMobj(FD,W);

S = klargestconncomps(S,3);
S = trunk(S);
CS  = STREAMobj2cell(S);

z1 = getnal(CS{1},DEM);
z2 = getnal(CS{2},DEM);
z3 = getnal(CS{3},DEM);

% ## plot stream networks
h1 = subplot(2,1,1);
hold(h1,'on');  
axes(h1);
imageschs(DEM);
hold on
clrs=jet(3);
% plot0(S,'k')
for r=1:numel(CS)
axes(h1)
% plot(CS{r},'color',clrs(r,:));
% axes(h2)
plot0(CS{r},'color',clrs(r,:));
end
hold off

% ## plot river profile
h2 = subplot(2,1,2);
hold(h2,'on');   
for r=1:numel(CS)
axes(h2)
% plot(CS{r},'color',clrs(r,:));
% axes(h2)
plotdz(CS{r},DEM,'color',clrs(r,:));
end
hold(h2,'off');
filename = "Stream_networks";
savefig(filename);

%figure();
%plot(CS{3}.distance,z3);


% ## plot the specific river profile
for r=1:numel(CS)
% plot figure
figure();    
filename_fig = "river"+ num2str(r);
% plot(CS{r},'color',clrs(r,:));
% axes(h2);
plotdz(CS{r},DEM,'color',clrs(r,:));
savefig(filename_fig);

% save data
h_line=get(gca,'Children');%get linehandles
distance=get(h_line,'Xdata');
elevation=get(h_line,'Ydata');
% figure();
% plot(xdata,ydata);
filename_data = "river"+ num2str(r) + '.mat';
save(['Data_rivers/', filename_data],'distance','elevation')
end


