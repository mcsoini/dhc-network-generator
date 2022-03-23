addpath(genpath('C:\Users\ashreeta\Documents\Martin\OSM\'))

stind = 2
cr_ind = 3

heat_dmnd = 131187*1e3 % GJ

%%
bnd_stt = m_shaperead('bundeslaender_borders')

state_name = bnd_stt.NAME

cr_states = { ...
'baden-wuerttemberg', ...
'bayern', ...
'berlin', ...
'brandenburg', ...
'bremen', ...
'hamburg', ...
'hessen', ...
'mecklenburg-vorpommern', ...
'niedersachsen', ...
'nordrhein-westfalen', ...
'rheinland-pfalz', ...
'saarland', ...
'sachsen', ...
'sachsen-anhalt', ...
'schleswig-holstein', ...
'thueringen', ...
};

state_border_index=[10, 16, 6, 8, 12, 15, 11, 1, 14, 13, 3, 5, 7, 9, 2, 4];

% for istate = 1:length(bnd_stt.NAME)
%     bnd_stt.NAME{istate}
%     
%     bnd_stt_arr = [[nan nan]; bnd_stt.ncst{istate}(:,1:2); [nan nan]]
%     
%     plot(bnd_cr(:,1),bnd_cr(:,2))
%     hold on 
%     plot(bnd_stt_arr(:,1),bnd_stt_arr(:,2),'r')
%     pause
% end

%% WHOLE COUNTRY
BRD = m_shaperead('OSM_COUNTRY_BORDERS')

cntind = strfind(BRD.NAME,'Germany');
cntind = find(~cellfun(@isempty,cntind));

bnd_cr = BRD.ncst{cntind};
plot(bnd_cr(:,1),bnd_cr(:,2))


%% Import population data

grump_0 = importdata('euds00ag.asc',' ',6);
grump = grump_0;
popdens_0 = flipud(grump.data);
popdens = popdens_0;
popdens(popdens <=0) = 0;

grump_ncols = strsplit(grump.textdata{1});      grump_ncols = str2double(grump_ncols{2});
grump_nrows = strsplit(grump.textdata{2});      grump_nrows = str2double(grump_nrows{2});
grump_xllcorner = strsplit(grump.textdata{3});  grump_xllcorner = str2double(grump_xllcorner{2});
grump_yllcorner = strsplit(grump.textdata{4});  grump_yllcorner = str2double(grump_yllcorner{2});
grump_cellsize = strsplit(grump.textdata{5});   grump_cellsize = str2double(grump_cellsize{2});

% xvec = (0:+1:grump_ncols - 1) * grump_cellsize + grump_xllcorner + 0.5 * grump_cellsize;
% yvec = (grump_nrows - 1:-1:0) * grump_cellsize + grump_yllcorner + 0.5 * grump_cellsize;
xvec = (0:+1:grump_ncols - 1) * grump_cellsize + grump_xllcorner + 0.5 * grump_cellsize;
yvec = (0:+1:grump_nrows - 1) * grump_cellsize + grump_yllcorner + 0.5 * grump_cellsize;

[X,Y] = meshgrid(xvec,yvec.');

area = grump_cellsize^2 * 111 * 111.32*cos(Y*pi/180);

popnum = popdens.*area;

country_borders=[bnd_cr;[NaN,NaN]];

% Filter population data by census region
in=inpoly([X(1:end);Y(1:end)], [country_borders(:,1).';country_borders(:,2).']);

popdens(not(in)) = 0;
popnum(not(in)) = 0;

heatdens = popdens/max(popdens(:))*1000;
heatdens_shr = popdens/sum(popdens(:))*1000;
heatshare = popnum/sum(popnum(:))*1000000;

pcolorjw(xvec,yvec,log(heatshare))
hold on
for istate = 1:length(bnd_stt.NAME)
    bnd_stt.NAME{istate}
    
    bnd_stt_arr = [[nan nan]; bnd_stt.ncst{istate}(:,1:2); [nan nan]];
    
    plot(bnd_cr(:,1),bnd_cr(:,2))
    hold on 
    plot(bnd_stt_arr(:,1),bnd_stt_arr(:,2),'r')
    
end

% %% Calculate built area map
% [btd_mpd] = resample_globeland(xvec,yvec,popnum);
% 
% pcolorjw(xvec,yvec,log(btd_mpd))
% hold on
% pcolorjw(xvec,yvec,log(btd_mpd))
% plot(country_borders(:,1),country_borders(:,2))
% 
% % Calculate population density bins
% [popbinnum,popbins] = hist(popdens(popdens>0),1:10:ceil(max(popdens(:)/1000))*1000);
% 
% popunq = unique(popdens(:));

% %% Refine population count data
% popnum_ref = zeros(size(popnum));
% dbin = popbins(2) - popbins(1);
% for ibin = 1:length(popbins)
%     if mod(ibin,100) == 0
%         ibin
%     end
%     bin_ind = popdens >= popbins(ibin)-dbin/2 & popdens < popbins(ibin)+dbin/2;
%     btd_sum = sum(btd_mpd(bin_ind));
%     pop_sum = sum(popnum(bin_ind));
%     
%     popnum_ref(bin_ind) = pop_sum * btd_mpd(bin_ind)/btd_sum;
% end
% 
% popdens_ref = popnum_ref./area;
% popdens_ref(isnan(popdens_ref)) = 0;


%% Define population density bins

[heatbinnum,heatbins] = hist(heatdens(heatdens>0),0:1:ceil(max(heatdens(:)/1))*1);
dbin=heatbins(2)-heatbins(1);


for ibin = 1:length(heatbins)
    bincand{ibin} = find(heatdens>heatbins(ibin)-dbin/2 & heatdens<=heatbins(ibin)+dbin/2 & heatdens > 0);
    
    heat_share_bins(ibin) = sum(heatshare(bincand{ibin}));
end

%% 

formatSpec = '%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%.8f%s';

filename = 'GERMANY_output_1915_160511b'

T = readtable([filename '.dat'],'Delimiter',',', ...
    'Format',formatSpec)

for ipx = 1:length(T.bin_nr)
    ipx
    A = round(T.heatdens(ipx)*1e8)
    B = round(heatdens(T.geo_index(ipx))*1e8)
    if round(T.heatdens(ipx)*1e8 )== round(heatdens(T.geo_index(ipx))*1e8)
        heatdens_shr_list(ipx) = heatdens_shr(T.geo_index(ipx));
        heatshare_px(ipx) = heatshare(T.geo_index(ipx));
    else
        'error'
        [T.heatdens(ipx)*1e8 , round(heatdens(T.geo_index(ipx))*1e8)]
        pause
    end
end

heatshare_px_table = array2table(heatshare_px.', 'VariableNames',{'heatshare_pixel'})

T2 = [T,heatshare_px_table]

writetable(T2, [filename, '_new.dat'],'Delimiter',',')


%% 
timstr=datestr(datetime('now'),'_HHMM_yymmdd');
out_file_name = ['GERMANY_output',timstr,'.dat'];
out_file_name = 'GERMANY_output_1915_160511b.dat'

%coord_blacklist = [3018296, 3016424,3014553]
coord_blacklist = []

min_length_cr = cell(1,length(heatbins));
con_length_cr = cell(1,length(heatbins));
cpcn_length_cr = cell(1,length(heatbins));
NSTRV_cr = cell(1,length(heatbins));
NBGSV_cr = cell(1,length(heatbins));
NSTR_cr = cell(1,length(heatbins));
NBGS_cr = cell(1,length(heatbins));
XCRD_cr = cell(1,length(heatbins));
YCRD_cr = cell(1,length(heatbins));
HEAT_cr = cell(1,length(heatbins));
STAT_cr = cell(1,length(heatbins));
ERRO_cr = cell(1,length(heatbins));

%out_file_name = 'output_1918_160331h.dat'

% fileID = fopen(out_file_name,'w');
% formatSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n';
% csv_clms = {'bin_nr', 'pixel_in_bin_nr', 'pixel_valid_nr','min_length_cr','con_length_cr','cpcn_length_cr','NSTRV_cr','NBGSV_cr','NSTR_cr','NBGS_cr','XCRD_cr','YCRD_cr','North','South','West','East','heat_share_bins','heatdens','STAT_cr','geo_index','pixel_area','ERRO_cr'}
% fprintf(fileID,formatSpec,csv_clms{:});
% fclose(fileID)

rng(123123)

for istate = 7:length(cr_states)
    BGS = m_shaperead([cr_states{istate}, '_buildings'])
    RDS = m_shaperead([cr_states{istate}, '_roads'])

    bnd_stt_index = state_border_index(istate);

    bnd_stt_arr = [[nan nan]; bnd_stt.ncst{bnd_stt_index}(:,1:2); [nan nan]];
    
    % find pixels within the state
    in=inpoly([X(1:end);Y(1:end)], [bnd_stt_arr(:,1).';bnd_stt_arr(:,2).']);

    heatshare_stt = heatshare;
    heatshare_stt(not(in)) = 0;
    heatdens_stt = heatdens;
    heatdens_stt(not(in)) = 0;

    pcolorjw(xvec,yvec,log(heatshare_stt))
    hold on 
    plot(bnd_stt.ncst{bnd_stt_index}(:,1:2),bnd_stt.ncst{bnd_stt_index}(:,2))

    % get state pixels in heatbins
    bincand_stt=cell(length(heatbins),1);
    for ibin = 1:length(heatbins)
        bincand_stt{ibin} = (find(heatdens_stt>heatbins(ibin)-dbin/2 & heatdens_stt<heatbins(ibin)+dbin/2  & heatdens_stt > 0));
        [ibin length(bincand_stt{ibin})];
        
        randind=randperm(length(bincand_stt{ibin}),length(bincand_stt{ibin}));
        binslct_stt{ibin} = bincand_stt{ibin}(randind);
    end
    
    %%
    time_dur = cell(1,length(binslct_stt));
    min_length = cell(1,length(binslct_stt));
    con_length = cell(1,length(binslct_stt));
    cpcn_length = cell(1,length(binslct_stt));
    NSTRV = cell(1,length(binslct_stt));
    NBGSV = cell(1,length(binslct_stt));
    NSTR = cell(1,length(binslct_stt));
    NBGS = cell(1,length(binslct_stt));
    XCRD = cell(1,length(binslct_stt));
    YCRD = cell(1,length(binslct_stt));
    HEAT = cell(1,length(binslct_stt));
    ERRO = cell(1,length(binslct_stt));

    frac = 0;
    
    if istate == 7
        idns_0 = 474
    else
        idns_0 = 1
    end
    
    for idns = idns_0:length(binslct_stt)

        indslct = binslct_stt{idns};
        indcand = bincand_stt{idns};

        frac = frac + sum(heatshare_stt(indcand))/sum(sum(heatshare_stt));

        iind = 1;
        cntval = 1;
        while cntval <= 1000000 & iind <= length(bincand_stt{idns}) & not(ismember(binslct_stt{idns}(iind),coord_blacklist))

    %         figure(543)
    %         pcolorjw(X,Y,log(popdata))
    %         hold on
    %         plot(X(indcand), Y(indcand),'.r')
    %         plot(X(indcand), Y(indcand),'.r')
    %         plot(X(binslct{idns}(iind)), Y(binslct{idns}(iind)),'og')
    %         hold off
    %         drawnow

            [istate iind length(binslct_stt{idns}) cntval idns length(binslct_stt)]

            coord = [...
                Y(binslct_stt{idns}(iind)) + grump_cellsize * 0.5;
                X(binslct_stt{idns}(iind)) + grump_cellsize * 0.5;
                Y(binslct_stt{idns}(iind)) - grump_cellsize * 0.5;
                X(binslct_stt{idns}(iind)) - grump_cellsize * 0.5];

            deg2metx=111;
            deg2mety=111.32*cos(Y(binslct_stt{idns}(iind))*pi/180);

            XCRD{idns}(cntval) = X(binslct_stt{idns}(iind));
            YCRD{idns}(cntval) = Y(binslct_stt{idns}(iind));
            HEAT{idns}(cntval) = heatshare_stt(binslct_stt{idns}(iind));
            STAT{idns}(cntval) = istate;
            ERRO{idns}{cntval} = [];

            try
                [time_dur{idns}(cntval), min_length{idns}(cntval), con_length{idns}(cntval), cpcn_length{idns}(cntval), NSTRV{idns}(cntval), NBGSV{idns}(cntval), NSTR{idns}(cntval), NBGS{idns}(cntval)] = create_network(coord, RDS, BGS, deg2mety,deg2metx);
                ERRO{idns}{cntval} = 'no error';
            catch ME
                time_dur{idns}(cntval) = NaN;
                min_length{idns}(cntval) = NaN;
                con_length{idns}(cntval) = NaN;
                cpcn_length{idns}(cntval) = NaN;
                NSTRV{idns}(cntval) = NaN;
                NBGSV{idns}(cntval) = NaN;
                NSTR{idns}(cntval) = NaN;
                NBGS{idns}(cntval) = NaN;
                ERRO{idns}{cntval} = strrep([ME.message,''],',','');
            end

            if not(strcmp(ERRO{idns}{cntval},'no streets no buildings: break')) & ...
               not(strcmp(ERRO{idns}{cntval},'no streets: break')) & ...
               not(strcmp(ERRO{idns}{cntval},'no buildings: break'))
                fileID = fopen(out_file_name,'a');
                formatSpec = '%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%d,%d,%.8f,%s\n';
                csv_clms = {idns, iind, cntval, min_length{idns}(cntval),con_length{idns}(cntval),cpcn_length{idns}(cntval),NSTRV{idns}(cntval),NBGSV{idns}(cntval),NSTR{idns}(cntval),NBGS{idns}(cntval),XCRD{idns}(cntval),YCRD{idns}(cntval),coord(1),coord(3),coord(4),coord(2),heat_share_bins(idns),heatdens_stt(binslct_stt{idns}(iind)),STAT{idns}(cntval),binslct_stt{idns}(iind),area(binslct_stt{idns}(iind)),ERRO{idns}{cntval}}
                fprintf(fileID,formatSpec,csv_clms{:});
                fclose(fileID)
            end

            if not(isnan(min_length{idns}(cntval)))
                cntval = cntval + 1;
            end

            iind = iind +1 ;        
        end
    end
    min_length_vec = []
    con_length_vec = []
    cpcn_length_vec = []
    NSTRV_vec = []
    NBGSV_vec = []
    NSTR_vec = []
    NBGS_vec = []
    HEAT_vec = []
    for i=1:length(min_length)
        min_length_vec = [min_length_vec; min_length{i}.'];
        con_length_vec = [con_length_vec; con_length{i}.'];
        cpcn_length_vec = [cpcn_length_vec; cpcn_length{i}.'];
        NSTRV_vec = [NSTRV_vec; NSTRV{i}.'];
        NBGSV_vec = [NBGSV_vec; NBGSV{i}.'];
        NSTR_vec = [NSTR_vec; NSTR{i}.'];
        NBGS_vec = [NBGS_vec; NBGS{i}.'];
        HEAT_vec = [HEAT_vec; HEAT{i}.'];

    end

    for ibin = 1:length(heatbins)
        min_length_cr{ibin} = [min_length_cr{ibin} min_length{ibin} ];
        con_length_cr{ibin}=[con_length_cr{ibin},con_length{ibin}];
        cpcn_length_cr{ibin}=[cpcn_length_cr{ibin},cpcn_length{ibin}];
        NSTRV_cr{ibin}=[NSTRV_cr{ibin},NSTRV{ibin}];
        NBGSV_cr{ibin}=[NBGSV_cr{ibin},NBGSV{ibin}];
        NSTR_cr{ibin}=[NSTR_cr{ibin},NSTR{ibin}];
        NBGS_cr{ibin}=[NBGS_cr{ibin},NBGS{ibin}];
        XCRD_cr{ibin}=[XCRD_cr{ibin},XCRD{ibin}];
        YCRD_cr{ibin}=[YCRD_cr{ibin},YCRD{ibin}];
        HEAT_cr{ibin}=[HEAT_cr{ibin},HEAT{ibin}];
        STAT_cr{ibin}=[STAT_cr{ibin},istate*ones(size(HEAT{ibin}))];
    end
    
    plot(HEAT_vec, (min_length_vec+con_length_vec+cpcn_length_vec)./HEAT_vec, 'o')
    drawnow
end
