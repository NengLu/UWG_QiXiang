DEM = GRIDobj('../../Data/dem/TRR_Earth2014.TBI2014.1min.order10800.tif');
%DEM = GRIDobj('../../Data/dem/TRR_Earth2014.BED2014.1min.order10800.tif');
DEM.refmat(3,1) = 88.0083333;
DEM.refmat(3,2) = 35.9916667;

%DEM = GRIDobj('../../Data/dem/TRR_SRTM.3s.size1200.tif');
%DEM.refmat(3,1) = 87.9995833;
%DEM.refmat(3,2) = 36.0004167;

th_discharge = 500; 
DEMf = fillsinks(DEM);
FD = FLOWobj(DEMf);

A = flowacc(FD);
W = A>th_discharge;

S = STREAMobj(FD,W);
S = klargestconncomps(S,5);
S = trunk(S);
CS  = STREAMobj2cell(S);

z1 = getnal(CS{1},DEM);
z2 = getnal(CS{2},DEM);
z3 = getnal(CS{3},DEM);

h1 = subplot(2,1,1);
hold(h1,'on');  
axes(h1);
imageschs(DEM);
hold on
clrs=jet(numel(CS));
%plot0(S,'k')
for r=1:numel(CS)
axes(h1)
%plot(CS{r},'color',clrs(r,:));
% axes(h2)
plot0(CS{r},'color',clrs(r,:));
end
hold off

h2 = subplot(2,1,2);
hold(h2,'on');   
for r=1:numel(CS)
axes(h2)
plot0(CS{r},'color',clrs(r,:));
%axes(h2)
%plotdz(CS{r},DEM,'color',clrs(r,:));

river_x = CS{r}.x;
river_y = CS{r}.y;
river_d = CS{r}.distance;
fname = "river_xyd"+num2str(r)+'.mat';
save(fname,'river_x','river_y','river_d');
end
hold(h2,'off')

%figure();
%plot(CS{3}.distance,z3);