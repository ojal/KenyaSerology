%% Plot a spatial map of peak timing

S = shaperead('ken_admbnda_adm1_iebc_20191031.shp');
[~,index] = sortrows({S.ADM1_EN}.'); S = S(index); clear index;
for i = 1:47
      S(i).ID = i;
end
load('plotsforpaper/fittedpeaktimesbycounty.mat');
load('plotsforpaper/datainferredpeaktimesbycounty.mat');
%% Combine the directly fitted peak times with the

peaktimes(peaktimes_fitted > 0 ) = peaktimes_fitted(peaktimes_fitted > 0 );
peaktimes(peaktimes_data_inferred> 0 ) = peaktimes_data_inferred(peaktimes_data_inferred > 0 );


%%
% min_peak_time = min(peaktimes(peaktimes > 0) );
%max_peak_time = max(peaktimes(peaktimes > 0) );

min_peak_time = 70; %aligns with May 1st
max_peak_time = 254; %aligns with Nov 1st

duration = max_peak_time - min_peak_time;

%%
n = round(max_peak_time - min_peak_time + 1);
map =parula(n);
% map = flipud(cool(n));


%%
greyshade = [0.7 0.7 0.7];
if peaktimes(1) < 0
    ColoredConstituencies = makesymbolspec('Polygon',{ 'ID',1,'FaceColor',greyshade} );
else
    I = ceil(peaktimes(1) - min_peak_time);
    ColoredConstituencies = makesymbolspec('Polygon',{ 'ID',1,'FaceColor',map(I,:)} );
end
    
    
 for i = 2:47
        if peaktimes(i) < 0
            ColoredConstituencies.FaceColor{i,1} = 'ID';
            ColoredConstituencies.FaceColor{i,2} = i;
            ColoredConstituencies.FaceColor{i,3} = greyshade;
        else
               I = ceil(peaktimes(i) - min_peak_time +0.1);
               ColoredConstituencies.FaceColor{i,1} = 'ID';
                ColoredConstituencies.FaceColor{i,2} = i;
                ColoredConstituencies.FaceColor{i,3} = map(I,:);
        end
 end
 
 for i=1:47
    if peaktimes_fitted(i) < 0
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
% mapshow(S,'LineStyle','none','SymbolSpec',ColoredConstituencies);
h = mapshow(S,'SymbolSpec',ColoredConstituencies);

for i=1:length(h.Children)
    if ToBeHashed(i)
        candystripe(h.Children(i),'Width',1); drawnow;
    end
end

cb = colorbar;
cb.Location = 'westoutside';
cb.Ticks = ([datenum(2020,5,1),datenum(2020,6,1),datenum(2020,7,1),datenum(2020,8,1),datenum(2020,9,1),datenum(2020,10,1),datenum(2020,11,1)] - datenum(2020,5,1))/duration;

cb.TickLabels = {'May','June','July','Aug','Sept','Oct','Nov'} ;
cb.Label.String = 'Date of peak incidence';
cb.FontSize = 28;



axes(h_Nairobi);
    mapshow(S(30),'SymbolSpec',ColoredConstituencies);
    title('Nairobi')
    
   axes(h_mombasa);
    mapshow(S(28),'SymbolSpec',ColoredConstituencies);
    title('Mombasa')  

%%

%  dim = [.496 .227 .025 .025];
annotation('ellipse',[0.49365625,0.227008149010477,0.0211875,0.026155995343423]);
annotation('ellipse',[0.356345911949686,0.415263748597082,0.032018867924528,0.032547699214368]);


%%

annotation('arrow',[0.5140625,0.73046875],[0.238813736903376,0.2782305005]);
annotation('arrow',[0.38671875,0.7203125],[0.437882421420256,0.707799767171129]);  

