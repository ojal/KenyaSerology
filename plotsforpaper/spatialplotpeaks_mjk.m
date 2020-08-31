%Plot a spatial map of peak timing

S = shaperead('County.shp');
%[~,index] = sortrows({S.COUNTY}.'); S = S(index); clear index;
for i = 1:47
      S(i).ID = i;
end
load('plotsforpaper/peaktimesbycounty.mat');
%%
map =(parula(100));
%%
I = ceil(rand()*100);
ColoredConstituencies = makesymbolspec('Polygon',{ 'ID',1,'FaceColor',map(I,:)} );
 for i = 2:47
        I = ceil(rand()*100);
        ColoredConstituencies.FaceColor{i,1} = 'ID';
        ColoredConstituencies.FaceColor{i,2} = i;
        ColoredConstituencies.FaceColor{i,3} = map(I,:);
 end
 
for i=2:47
    if rand()<0.3
        ToBeHashed(i)=1;
    else
        ToBeHashed(i)=0;
    end
end
    
 %%
colormap(map);
figure(1);
clf;
%Create three axes
h_kenya = axes;
h_kenya.Position =  [0.100 0.100 0.55 0.85];
h_kenya.FontSize = 18;
h_kenya.XAxis.Visible = 'off';
h_kenya.YAxis.Visible = 'off';


h_Nairobi = axes;
h_Nairobi.Position =  [0.7 0.600 0.25 0.35];
h_Nairobi.FontSize = 18;
h_Nairobi.XAxis.Visible = 'off';
h_Nairobi.YAxis.Visible = 'off';

h_mombasa = axes;
h_mombasa.Position =  [0.7 0.100 0.25 0.35];
h_mombasa.FontSize = 18;
h_mombasa.XAxis.Visible = 'off';
h_mombasa.YAxis.Visible = 'off';

axes(h_kenya);
h=mapshow(S,'LineStyle','none','SymbolSpec',ColoredConstituencies);

% drawnow;

for i=1:length(h.Children)
    if ToBeHashed(i)
        candystripe(h.Children(i),'Width',1); drawnow;
    end
end

cb = colorbar;
cb.Location = 'westoutside';
cb.Ticks = 0:0.5:1;
cb.TickLabels = {'1st June','1st July','1st Aug'} ;
cb.Label.String = 'Peak date';
cb.FontSize = 18;
%     cb.Ticks = [log10(1)/log10(max_inc);log10(11)/log10(max_inc);log10(101)/log10(max_inc);log10(1001)/log10(max_inc);log10(10001)/log10(max_inc)];
%     cb.TickLabels = [0;10;100;1000;10000];
%     cb.Limits = [0 1];

% 
% axes(h_Nairobi);
%     mapshow(S(30),'SymbolSpec',ColoredConstituencies);
%     title('Nairobi')
%     
%    axes(h_mombasa);
%     mapshow(S(28),'SymbolSpec',ColoredConstituencies);
%     title('Mombasa')  

    
%     axes(h_kenya);
%     cb = colorbar;
    