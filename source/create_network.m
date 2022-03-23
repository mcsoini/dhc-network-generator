function [time_dur, nw_length, con_length, cpcn_length, NSTRV, NBGSV, NSTR, NBGS] = create_network(coord, RDS, BGS, is_bgs, deg2mety,deg2metx)
tic
% 
% clearvars -except BGS RDS area  coord frac grump_cellsize grump_ncols grump_nrwos grump_xllcorner grump_yllcorner iind iunq max_ind popcnt popdata popdataunq RDS X xvec Y yvec
%% 
% coord
% RDS
% BGS
% deg2mety
% deg2metx
% % 


limN = coord(1) ;
limE = coord(2) ;
limS = coord(3) ;
limW = coord(4) ;

%%%%%%%%%%
% FILTER %

% STREETS

% Bridge
%str_isbrdg = cat(1,ismember({RDS.bridge{:}},'yes')).';

% Tunnel
% str_istnnl = cat(1,ismember({RDS.tunnel{:}},'yes')).';

% Proper street
%str_isstr = cat(1,(ismember({RDS.highway{:}},{'pedestrian','living_street','trunk','service','footway','pedestrian','primary','secondary','tertiary','unclassified','residential','trunk_link','primary_link','secondary_link','tertiary_link'}))).';

str_slct = ...
   ((RDS.mbr(:,1) > limW & RDS.mbr(:,1) < limE) | (RDS.mbr(:,3) < limE & RDS.mbr(:,3) > limW)) &...
   ((RDS.mbr(:,4) < limN & RDS.mbr(:,4) > limS) | (RDS.mbr(:,2) > limS & RDS.mbr(:,2) < limN)); 
%     str_istnnl == 0;
%     str_isstr == 1;

% BUILDINGS

% Is building?
% bgs_bdg = cat(1,strfind({BGS.building{:}},'yes'));
% bgs_bdg(cellfun(@isempty,bgs_bdg)) = {0};
% bgs_isbgs = cat(1,bgs_bdg{:});
% 
is_bgs_bool = zeros(size(BGS.mbr(:,1)));
is_bgs_bool(is_bgs) = 1;


bgs_slct = ...
    BGS.mbr(:,1) > limW & ...
    BGS.mbr(:,2) > limS & ...
    BGS.mbr(:,3) < limE & ...
    BGS.mbr(:,4) < limN & ...
    is_bgs_bool;


%%
str_slct_coord = {RDS.ncst{str_slct}};
bgs_slct_coord = {BGS.ncst{bgs_slct}};

% convert coordinates to meters and round to integer millimeters
bgs_slct_coord_mm = cell(size(bgs_slct_coord));
str_slct_coord_mm = cell(size(str_slct_coord));
for ibgs = 1:length(bgs_slct_coord)
     bgs_slct_coord_mm{ibgs}(:,1) = round((bgs_slct_coord{ibgs}(:,1) - limW) * deg2metx*1000 * 1000);
     bgs_slct_coord_mm{ibgs}(:,2) = round((bgs_slct_coord{ibgs}(:,2) - limS) * deg2mety*1000 * 1000);
end
for istr = 1:length(str_slct_coord_mm)
     str_slct_coord_mm{istr}(:,1) = round((str_slct_coord{istr}(:,1) - limW) * deg2metx*1000 * 1000);
     str_slct_coord_mm{istr}(:,2) = round((str_slct_coord{istr}(:,2) - limS) * deg2mety*1000 * 1000);
end
    
NSTRV = numel(vertcat(str_slct_coord_mm{:}));
NBGSV = numel(vertcat(bgs_slct_coord_mm{:}));

NSTR = numel(str_slct_coord_mm);
NBGS = numel(bgs_slct_coord_mm);


bgs_slct_coord_nan={};
str_slct_coord_nan={};
str_slct_coord_nan{max(1,size(str_slct_coord_mm,1)),max(1,size(str_slct_coord_mm,2))} = [];
bgs_slct_coord_nan{max(1,size(bgs_slct_coord_mm,1)),max(1,size(bgs_slct_coord_mm,2))} = [];



% Add NaN cells and round street and building coordinates
for istr_slct=1:length(str_slct_coord_mm)
    str_slct_coord_nan{1,istr_slct} = str_slct_coord_mm{istr_slct};
    str_slct_coord_nan{2,istr_slct} = [nan,nan];
end
for ibgs_slct=1:length(bgs_slct_coord_mm)
    bgs_slct_coord_nan{1,ibgs_slct} = bgs_slct_coord_mm{ibgs_slct};
    bgs_slct_coord_nan{2,ibgs_slct} = [nan,nan];
end

% Separated array of street and building coordinates
str_slct_coord_array_nan = [[NaN, NaN]; cat(1,str_slct_coord_nan{:})];
bgs_slct_coord_array_nan = [[NaN, NaN]; cat(1,bgs_slct_coord_nan{:})];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1231)
% plot(str_slct_coord_array_nan(:,1),str_slct_coord_array_nan(:,2),'k.-')
% hold on
% plot(bgs_slct_coord_array_nan(:,1),bgs_slct_coord_array_nan(:,2),'r-')
% rectangle('Position',[0 0 (limE-limW)*deg2metx*1e6 (limN-limS)*deg2mety*1e6])
% hold off
% drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(str_slct_coord_mm) == 0 || length(bgs_slct_coord_mm) == 0
    time_dur = NaN;
    nw_length = NaN;
    con_length = NaN;
    cpcn_length = NaN;
    NSTRV = NaN;
    NBGSV = NaN;
    NSTR = NaN;
    NBGS = NaN;
    if length(str_slct_coord_mm) == 0 && length(bgs_slct_coord_mm) == 0
        error('no streets, no buildings: break')
    elseif length(str_slct_coord_mm) < 2
        error('no streets: break')
    elseif length(bgs_slct_coord_mm) < 2
        error('no buildings: break')
    end
end



% Create list of unique original street vertices before modifying the
% streets
str_vertices_0 = unique(cat(1,str_slct_coord_mm{1,:}),'rows');

%%
str_vertices_cnnct = str_vertices_0;

nbr_str = zeros(1,length(bgs_slct_coord_mm));
x_str_connect = zeros(1,length(bgs_slct_coord_mm));
y_str_connect = zeros(1,length(bgs_slct_coord_mm));
x_bgs_connect = zeros(1,length(bgs_slct_coord_mm));
y_bgs_connect = zeros(1,length(bgs_slct_coord_mm));
str_min_d = zeros(1,length(bgs_slct_coord_mm));
ind_str_connect = zeros(1,length(bgs_slct_coord_mm));
ind_str_pos = zeros(1,length(bgs_slct_coord_mm));

cnnct_vert1 = zeros(length(bgs_slct_coord_mm),2);

cnnct_vert2 = zeros(length(bgs_slct_coord_mm),2);

length(bgs_slct_coord_mm);

cnt_zro=0;

% Loop over buildings to find shortest connection to road for each building
for ibgs = 1:length(bgs_slct_coord_mm)
     if mod(ibgs,100)==0
        [ibgs length(bgs_slct_coord_mm)]
     end

    % index of the street the building is connected to
    ind_str_connect(ibgs) = 0;
    % x-coordinate of the street connection point
    x_str_connect(ibgs) = -1;
    % y-coordinate of the street connection point
    y_str_connect(ibgs) = -1;
    % x-coordinate of the building connection point
    x_bgs_connect(ibgs) = -1;
    % y-coordinate of the building connection point
    y_bgs_connect(ibgs) = -1;
    % distance from the building to the connection point
    str_min_d(ibgs) = 10^9;
    % number of streets within the search radius
    nbr_str(ibgs) = 0;
    
    % radius threshold for the initial street search
    radius_thrshld = 60*1e3;

    while nbr_str(ibgs) == 0 % loop until one street is within the initial street search radius, otherwise increase radius_thrshld
 
        % loop over all streets
        for istr = 1:length(str_slct_coord_mm)
        
        %calculate the center point of the building
        bgs_slct_coord_cntr = sum(bgs_slct_coord_mm{1,ibgs},1)/length(bgs_slct_coord_mm{1,ibgs}(:,1));
        
            % check if any of the vertices of the current street is close
            % enough to the building center point
            if min(sum((str_slct_coord_mm{1,istr} - repmat(bgs_slct_coord_cntr,size(str_slct_coord_mm{1,istr},1),1)).^2,2)) < radius_thrshld^2

                % increase counter of streets within search radius
                nbr_str(ibgs) = nbr_str(ibgs) + 1 ;

                % assemble building and street polygon coordinates as an
                % input to the function min_dist_between_two_polygons
                VB.x = bgs_slct_coord_mm{1,ibgs}(:,1);
                VB.y = bgs_slct_coord_mm{1,ibgs}(:,2);
                ST.x = str_slct_coord_mm{1,istr}(:,1);
                ST.y = str_slct_coord_mm{1,istr}(:,2);

                % calculate minimal distance between building and street
                [STx, STy, BDx, BDy, min_d] = min_dist_between_two_polygons(VB,ST);
%                  figure(141235)
%                 plot(str_slct_coord_mm{1,istr}(:,1),str_slct_coord_mm{1,istr}(:,2),'ko-')
%                 hold on 
%                 plot([STx BDx],[ STy BDy],'rx-')
%                 plot(bgs_slct_coord_mm{1,ibgs}(:,1),bgs_slct_coord_mm{1,ibgs}(:,2),'b.-')

%                 drawnow
%                 pause 
                
                % check if the current street is the closest so far
                if min_d < str_min_d(ibgs)
                    str_min_d(ibgs) = min_d;        % list with minimum connection distances
                    ind_str_connect(ibgs) = istr;   % list with streets the buildings are connected to
                    x_str_connect(ibgs) = STx;      % list with coordinates of the street connection points
                    y_str_connect(ibgs) = STy;      % list with coordinates of the street connection points
                    x_bgs_connect(ibgs) = BDx;      % list with coordinates of the building connection points
                    y_bgs_connect(ibgs) = BDy;      % list with coordinates of the building connection points
                end
            end
        end
        radius_thrshld = radius_thrshld * 2; % double the street search radius if no street was close enough
    end

%     plot(VB.x,VB.y,'.b-')
%     hold on
%     plot(str_slct_coord_mm{1,ind_str_connect(ibgs)}(:,1),str_slct_coord_mm{1,ind_str_connect(ibgs)}(:,2),'o-r')
%     plot([x_str_connect(ibgs),x_bgs_connect(ibgs)],[y_str_connect(ibgs),y_bgs_connect(ibgs)],'x-g')
%     rectangle('Position',[0 0 (limE-limW)*deg2metx*1e6 (limN-limS)*deg2mety*1e6])
%     text(0.5*(x_str_connect(ibgs)+x_bgs_connect(ibgs)), 0.5*(y_str_connect(ibgs)+y_bgs_connect(ibgs)) ,num2str(str_min_d(ibgs)/1000,2))
% 
%     drawnow
%     pause


    % Check position of the street connection point with respect to the
    % street vertices


    
    vert_street = str_slct_coord_mm{1,ind_str_connect(ibgs)};

    r_mix_x = (x_str_connect(ibgs) - vert_street(:,1))./([vert_street(2:end,1);0] - vert_street(:,1));
    r_mix_y = (y_str_connect(ibgs) - vert_street(:,2))./([vert_street(2:end,2);0] - vert_street(:,2));
%     [r_mix_x r_mix_y]
    
%     plot(vert_street(:,1),vert_street(:,2),'o-')
%     hold on
%     plot(x_str_connect(ibgs), y_str_connect(ibgs), 'x')
%     drawnow 
%     pause
%     hold off
    
    % if the coordinates of the two street vertices are equal r_mix_x
    % and/or r_mix_y are NaN
    if sum(isnan(r_mix_x)) == 0 & sum(isnan(r_mix_y)) == 0
         % vertices with both x and y coordinate equal the connection point
        [~,vert_match] = min(abs(r_mix_x - r_mix_y));
    elseif sum(isnan(r_mix_x)) > 0 & sum(isnan(r_mix_y)) > 0 | ...
            sum(isinf(r_mix_x)) > 0 & sum(isnan(r_mix_y)) > 0 | ...
            sum(isnan(r_mix_x)) > 0 & sum(isinf(r_mix_y)) > 0
        vert_match = find(isnan(r_mix_x + r_mix_y));
        vert_match = vert_match(1);
    elseif sum(isnan(r_mix_x)) == 0 & sum(isnan(r_mix_y)) > 0
        vert_match = find(isnan(r_mix_y) & r_mix_x<=1 & r_mix_x>=0);
        vert_match = vert_match(1);
    elseif sum(isnan(r_mix_x)) > 0 & sum(isnan(r_mix_y)) == 0
        vert_match = find(isnan(r_mix_x) & r_mix_y<=1 & r_mix_y>=0);
        vert_match = vert_match(1);
    end
    
    % If connection point equals an existing vertex
    if (r_mix_x(vert_match) == 0 &  r_mix_y(vert_match) == 0) | (r_mix_x(vert_match) == 1 &  r_mix_y(vert_match) == 1)
        
        % get index of matching vertex
        vert_match_ind = (vert_match);
        %restrict to first if several identical vertices
        vert_match_ind = vert_match_ind(1);
        
        % set the first neighboring vertex equal the identical vertex
        c_vert1 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(vert_match_ind,:);
        % if the matching vertex is the first of the street...
        if vert_match_ind == 1
            % ...set the second neighbor equal the second street vertex
            c_vert2 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(vert_match_ind + 1, :);        
        elseif vert_match_ind == length(str_slct_coord_mm{1,ind_str_connect(ibgs)})
            c_vert2 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(vert_match_ind - 1, :)
        else
            % ...set the second neighbor equal the previous vertex
            c_vert2 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(vert_match_ind + 1, :);
            c_vert1 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(vert_match_ind, :);
        end
        
        ind_str_pos(ibgs) = vert_match_ind;

    % if the street has only two vertices
    elseif size(str_slct_coord_mm{1,ind_str_connect(ibgs)},1) == 2
        c_vert1 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(1,:);
        c_vert2 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(2,:);

        ind_str_pos(ibgs) = 1; % street vertex index before the connection point
    % if the street has more vertices and the connection point is in between    
    else
        
        
        ind_c_vert1 = vert_match;
        ind_c_vert2 = ind_c_vert1 + 1;

        c_vert1 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(ind_c_vert1,:);
        c_vert2 = str_slct_coord_mm{1,ind_str_connect(ibgs)}(ind_c_vert2,:);

        ind_str_pos(ibgs) = ind_c_vert1; % street vertex index before the connection point
    end

    %First neighboring vertex
    cnnct_vert1(ibgs,:) = c_vert1; 
    %Second neighboring vertex
    cnnct_vert2(ibgs,:) = c_vert2; 

    % make street connection point a new street coordinate
    % special cases here
    
    
%     'index'
%     [ind_str_pos(ibgs) size(str_slct_coord{1,ind_str_connect(ibgs)},1)]; 
%     
%     'first'
%     str_slct_coord{1,ind_str_connect(ibgs)}(1:ind_str_pos(ibgs),:); 
%     
%     'then'
%     [x_str_connect(ibgs) y_str_connect(ibgs)]; 
%     
%     'last'
%     str_slct_coord{1,ind_str_connect(ibgs)}(ind_str_pos(ibgs) + 1:end,:); 
    
    
    str_slct_coord_old = str_slct_coord_mm;

%     % else
%     if ind_str_connect(ibgs) == 9
%         plot(str_slct_coord{1,ind_str_connect(ibgs)}(:,1),str_slct_coord{1,ind_str_connect(ibgs)}(:,2),'o-b')
%         hold on
%     end
    
    str_slct_coord_mm{1,ind_str_connect(ibgs)} = ...
        [str_slct_coord_mm{1,ind_str_connect(ibgs)}(1:ind_str_pos(ibgs),:); ...
        [x_str_connect(ibgs) y_str_connect(ibgs)]; ...
        str_slct_coord_mm{1,ind_str_connect(ibgs)}(ind_str_pos(ibgs) + 1:end,:)];
    
%     if ind_str_connect(ibgs) == 9
%         plot(str_slct_coord{1,ind_str_connect(ibgs)}(:,1),str_slct_coord{1,ind_str_connect(ibgs)}(:,2),'x-r')
%         hold off
%         
%         str_slct_coord_old{1,ind_str_connect(ibgs)}
%         str_slct_coord{1,ind_str_connect(ibgs)}
%         
%         pause
%     end
%         
%     
 
    %x-coordinate of the new vertex
    %x_str_connect(ibgs)
    %y-coordinate of the new vertex
    %y_str_connect(ibgs)
    %street index of the new vertex
    %ind_str_connect(ibgs)

    
%     plot(str_slct_coord{1,ind_str_connect(ibgs)}(:,1),str_slct_coord{1,ind_str_connect(ibgs)}(:,2)*100000000,'o-')
%     hold on
%     plot(x_str_connect(ibgs), y_str_connect(ibgs), 'or')
%     plot(c_vert1(1),c_vert1(2),'xg')
%     plot(c_vert2(1),c_vert2(2),'xk')
%     hold off
% %     
%      drawnow

%     spy(str_adjmat)
%     drawnow
end


%%
% plot(str_slct_coord_array_nan(:,1)/100000000,str_slct_coord_array_nan(:,2)/100000000,'kx-')
% hold on
% 
% 
% for ibgs = 1:length(bgs_slct_coord)
%     ind_str_connect(ibgs)
%     if ind_str_connect(ibgs) ==264
%       %plot(str_slct_coord{1,ind_str_connect(ibgs)}(:,1),str_slct_coord{1,ind_str_connect(ibgs)}(:,2),'o-')
%       hold on
%       plot(x_str_connect(ibgs), y_str_connect(ibgs), 'or')
%       %plot(cnnct_vert1(ibgs,1),cnnct_vert1(ibgs,2),'xg')
%      % plot(cnnct_vert2(ibgs,1),cnnct_vert2(ibgs,2),'xg')
%      % 
%       plot([x_bgs_connect(ibgs),x_str_connect(ibgs)],[y_bgs_connect(ibgs),y_str_connect(ibgs)],'o-')
%       
%       plot(bgs_slct_coord{1,ibgs}(:,1), bgs_slct_coord{1,ibgs}(:,2))
%       
%       
%       drawnow
%       pause
%     end
% end

%%

str_vertices_1 = unique(cat(1,str_slct_coord_mm{1,:}),'rows');

% Delete vertices outside the area to avoid symmetry issues
str_vertices_1_outside = ...
    str_vertices_1(:,1) < 0 | ...
    str_vertices_1(:,2) < 0 | ...
    str_vertices_1(:,1) > (limE-limW)*deg2metx*1e6 | ...
    str_vertices_1(:,2) > (limN-limS)*deg2mety*1e6;
str_vertices_1(str_vertices_1_outside,:) = [];

str_slct_coord_nan_1=cell(size(str_slct_coord_mm,1),size(str_slct_coord_mm,2));

% Add NaN cells and round street and building coordinates
for istr_slct=1:length(str_slct_coord_mm)
    str_slct_coord_nan_1{1,istr_slct} = str_slct_coord_mm{istr_slct};
    str_slct_coord_nan_1{2,istr_slct} = [nan,nan];
end
% 
% plot(round(str_slct_coord{17}(:,1)*100000000),round(str_slct_coord{17}(:,2)*100000000),'.-')
% hold on
% plot(str_slct_coord_nan_1{1,17}(:,1),str_slct_coord_nan_1{1,17}(:,2),'.--')
% for istr = 1:length(str_slct_coord_nan_1)
%     plot(str_slct_coord_nan_1{1,istr}(:,1),str_slct_coord_nan_1{1,istr}(:,2),'.-')
%     text(sum(str_slct_coord_nan_1{1,istr}(:,1))/length(str_slct_coord_nan_1{1,istr}(:,1)),sum(str_slct_coord_nan_1{1,istr}(:,2))/length(str_slct_coord_nan_1{1,istr}(:,1)),num2str(istr,4))
% end

% Separated array of street and building coordinates
str_slct_coord_array_nan_1 = [[NaN, NaN]; cat(1,str_slct_coord_nan_1{:})];

% Located separators in the buildings list
str_slct_coord_array_isnan = sum(isnan(str_slct_coord_array_nan_1),2) ~= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       BUILD ADJACENCY MATRIX             %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate adjacency matrix of the street graph
str_adjmat = zeros(size(str_vertices_1,1));

for ivert = 1:size(str_vertices_1,1)
    if mod(ivert-1,100) == 0
        ['building adjacency matrix: ',num2str(ivert),'/',num2str(size(str_vertices_1,1))]
    end
    % search current vertex in the complete list of street vertices
    str_vertices_occ = (str_slct_coord_array_nan_1(:,1) == str_vertices_1(ivert,1) & str_slct_coord_array_nan_1(:,2) == str_vertices_1(ivert,2));
    
    % define adjacent vertex pair
    str_vertices_occ_prev = [str_vertices_occ(2:end);0];
    str_vertices_occ_next = [0; str_vertices_occ(1:end-1)];

    % combine adjacent vertex pair and make sure we stay within the street
    str_vertices_occ_adja_ind = (str_vertices_occ_prev + str_vertices_occ_next) == 1 & str_slct_coord_array_isnan == 0;
    
    % get adjacent vertex pair coordinates
    str_vertices_occ_adja = str_slct_coord_array_nan_1(str_vertices_occ_adja_ind,:);
    
    for iadj = 1:size(str_vertices_occ_adja,1)
        str_vertices_adja_dist = sqrt(sum((str_vertices_occ_adja(iadj,:) - str_vertices_1(ivert,:)).^2)) ;
        
        
        str_adjmat(str_vertices_1(:,1) == str_vertices_occ_adja(iadj,1) & str_vertices_1(:,2) == str_vertices_occ_adja(iadj,2), ivert) = str_vertices_adja_dist ;
      
%         figure(3241)
%         hold off
%         plot(str_slct_coord_array_nan(:,1),str_slct_coord_array_nan(:,2),'g')
%         hold on
%         plot(bgs_slct_coord_array_nan(:,1),bgs_slct_coord_array_nan(:,2))
%         drawnow
%         plot(str_vertices_1(:,1),str_vertices_1(:,2),'o')
%         hold on
%         for ivert = 1:length(str_vertices_1)
%             plotlist_neighb = find(str_adjmat(ivert,:));
% 
%             for ivert2 = 1:length(plotlist_neighb)
%                 if str_adjmat(ivert,plotlist_neighb(ivert2)) > 0
%                     plot([str_vertices_1(ivert,1),str_vertices_1(plotlist_neighb(ivert2),1)], [str_vertices_1(ivert,2),str_vertices_1(plotlist_neighb(ivert2),2)],'k-')
%                     drawnow
%                 end
%             end
%         end   
%         
%         if round(str_vertices_adja_dist*10000) == round(100459.646127189*10000)
%         'here'
%         
%             pause
%         
%         end

    end
end

cnt = 0;
for imat1 = 1:length(str_adjmat)
    for imat2 = 1:length(str_adjmat)
        if str_adjmat(imat1, imat2) ~=str_adjmat(imat2, imat1)
            cnt = cnt + 1
            [cnt imat1 imat2 str_adjmat(imat1, imat2) str_adjmat(imat2, imat1) ]
        end
    end
end
  

%%%%%%       BUILD ADJACENCY MATRIX             %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1231)
% hold on
% for i = 1:length(str_slct_coord_mm)
%     plot(str_slct_coord_mm{i}(:,1),str_slct_coord_mm{i}(:,2),'k.-')
% end
% plot(bgs_slct_coord_array_nan(:,1),bgs_slct_coord_array_nan(:,2))
% plot([str_vertices_1(39,1), str_vertices_1(41,1)],[str_vertices_1(39,2), str_vertices_1(41,2)],'o')
% rectangle('Position',[0 0 (limE-limW)*deg2metx*1e6 (limN-limS)*deg2mety*1e6])
% hold off
% drawnow


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       DELETE VERTICES W/O CONNECTION     %%%%%%%%%%%%%%%%%%%%%%%%%%%

str_vertices_d = ismember(round(str_vertices_1),round([x_str_connect.' y_str_connect.']));
str_vertices_d_log = str_vertices_d(:,1) == 0 & str_vertices_d(:,1) == 0;
str_vertices_d_ind = find(str_vertices_d_log);
str_vertices_d = str_vertices_1(str_vertices_d_log,:); % vertices to be deleted

% figure(13412)
% hold on

for ivdel1 = 1:length(str_vertices_d_ind)
    
    if mod(ivdel1-1,100) == 0
        ['deleting vertices: ',num2str(ivdel1),'/',num2str(length(str_vertices_d_ind))]
    end
    
    ind_neighb = find(str_adjmat(str_vertices_d_ind(ivdel1),:));
    
    % get all combinations of neighbors
    [ind_neighb_mesh_1,ind_neighb_mesh_2] = meshgrid(ind_neighb,ind_neighb);
    ind_neighb_comb = [ind_neighb_mesh_1(:) ind_neighb_mesh_2(:)];
    ind_neighb_comb = ind_neighb_comb(ind_neighb_comb(:,1) ~= ind_neighb_comb(:,2),:);
    
    if length(ind_neighb) > 1
        
        for icomb = 1:size(ind_neighb_comb,1)
            comb_dist = str_adjmat(str_vertices_d_ind(ivdel1), ind_neighb_comb(icomb,1)) + str_adjmat(str_vertices_d_ind(ivdel1), ind_neighb_comb(icomb,2));
            if str_adjmat(ind_neighb_comb(icomb,1),ind_neighb_comb(icomb,2)) == 0
             %   'new connection'
                
                str_adjmat(ind_neighb_comb(icomb,1),ind_neighb_comb(icomb,2)) = comb_dist;
            else
                %'minimizing existing distance'
                str_adjmat(ind_neighb_comb(icomb,1),ind_neighb_comb(icomb,2)) = min(str_adjmat(ind_neighb_comb(icomb,1),ind_neighb_comb(icomb,2)), comb_dist);
            end
            
        end
    end

    str_adjmat(str_vertices_d_ind(ivdel1),:) = zeros(1,length(str_adjmat));
    str_adjmat(:,str_vertices_d_ind(ivdel1)) = zeros(length(str_adjmat),1);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    if mod(ivdel1,1000)==0 || ivdel1 == length(str_vertices_d_ind)
%         figure(13412)
%         plot(str_vertices_d(ivdel1,1),str_vertices_d(ivdel1,2),'og')
%         hold off
%         drawnow

%pause

%         for ivert = 1:length(str_vertices_1)
%             plotlist_neighb = find(str_adjmat(ivert,:));
%             for ivert2 = 1:length(plotlist_neighb)
%                 if str_adjmat(ivert,plotlist_neighb(ivert2)) > 0
%                     plot([str_vertices_1(ivert,1),str_vertices_1(plotlist_neighb(ivert2),1)], [str_vertices_1(ivert,2),str_vertices_1(plotlist_neighb(ivert2),2)],'k-')
%                     hold on
% 
%                     % text(0.5*(str_vertices_1(ivert,1)+str_vertices_1(plotlist_neighb(ivert2),1)),0.5*(str_vertices_1(ivert,2)+str_vertices_1(plotlist_neighb(ivert2),2)),num2str(str_adjmat(ivert,plotlist_neighb(ivert2))/10^4,4))
%                 end
%             end
%         end
%         for ibgs = 1:length(x_bgs_connect)
%             plot([x_str_connect(ibgs)],[y_str_connect(ibgs)],'-ob')
%         end
%         plot(str_vertices_d(:,1),str_vertices_d(:,2),'xr')
%         drawnow
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
str_vertices_2 = str_vertices_1;
str_vertices_2(str_vertices_d_log,:) = [];

str_adjmat(not(str_vertices_d_log),not(str_vertices_d_log));

str_adjmat(str_vertices_d_log,:)=[];
str_adjmat(:,str_vertices_d_log)=[];

format longG
issymmetric(str_adjmat)
%%

size(str_adjmat)
str_adjmat_sparse = sparse(str_adjmat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(22343)
% plot(graph(str_adjmat))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_comp,ind_comp] = graphconncomp(str_adjmat_sparse,'Directed', false,'Weak', false);


comp_dist_vec = 0;

if n_comp >= 2
    nk_comp = nchoosek(1:n_comp,2);

    comp_dist=zeros(n_comp-1);
    if n_comp > 1
        for ink = 1:size(nk_comp,1)
            A=str_vertices_2(ind_comp==nk_comp(ink,1),:);
            B=str_vertices_2(ind_comp==nk_comp(ink,2),:);

            comp_dist(nk_comp(ink,1), nk_comp(ink,2)-1) = min(min(pdist2(A,B,'euclidean')));
        end
        comp_dist(comp_dist == 0) = 2*max(max(comp_dist));
        comp_dist_vec = min(comp_dist);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  GET CONNECTED COMPONENTS AND CALCULATE MINIMAL SPANNING TREE  %%%%%

str_adjmat_comp = cell(n_comp,1);

min_length = NaN*zeros(n_comp,1);
for icmp = 1:n_comp
    str_adjmat_comp{icmp} = str_adjmat_sparse(ind_comp==icmp,ind_comp==icmp);

    [str_adjmat_min,pred] = graphminspantree(str_adjmat_comp{icmp});
 %   view(biograph(str_adjmat_min,[],'ShowArrows','off','ShowWeights','on'))
 
    min_length(icmp) = full(sum(sum(str_adjmat_min)));
end
 
nw_length = sum(min_length);
con_length = sum(str_min_d);
cpcn_length = sum(comp_dist_vec);

%%%%%%  GET CONNECTED COMPONENTS AND CALCULATE MINIMAL SPANNING TREE  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% 
% % CHECK INTEGRITY OF ADJACENCY MATRIX
% figure(3241)
% plot(str_slct_coord_array_nan(:,1),str_slct_coord_array_nan(:,2),'g')
% hold on
% plot(bgs_slct_coord_array_nan(:,1),bgs_slct_coord_array_nan(:,2))
% drawnow
% plot(str_vertices_1(:,1),str_vertices_1(:,2),'o')
% hold on
% for ivert = 1:length(str_vertices_1)
%     plotlist_neighb = find(str_adjmat(ivert,:));
% 
%     for ivert2 = 1:length(plotlist_neighb)
%         if str_adjmat(ivert,plotlist_neighb(ivert2)) > 0
%             plot([str_vertices_1(ivert,1),str_vertices_1(plotlist_neighb(ivert2),1)], [str_vertices_1(ivert,2),str_vertices_1(plotlist_neighb(ivert2),2)],'k-')
%         end
%     end
% end
% % for ibgs = 1:length(x_bgs_connect)
% %     plot([x_str_connect(ibgs),x_bgs_connect(ibgs)]*100000000,[y_str_connect(ibgs),y_bgs_connect(ibgs)]*100000000,'-ob')
% % end
% % plot(str_vertices_d(:,1),str_vertices_d(:,2),'xr')

time_dur=toc
