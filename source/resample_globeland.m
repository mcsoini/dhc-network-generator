function [btd_mpd] = resample_globeland(xvec,yvec,popdata)
%%
UTMlist = [...
'10_30',
'10_35',
'10_40',
'10_45',
'11_30',
'11_35',
'11_40',
'11_45',
'12_30',
'12_35',
'12_40',
'12_45',
'13_30',
'13_35',
'13_40',
'13_45'];

btd_mpd = zeros(size(popdata));

for iutm = 1:length(UTMlist)
    BTD = ncread(['N',UTMlist(iutm,:),'_2010LC030_ATS_WGS.nc'],'Band1');

    BTD(BTD>0) = 1;
    BTD(BTD~=1) = 0;

    xbtd = ncread(['N',UTMlist(iutm,:),'_2010LC030_ATS_WGS.nc'],'lon');
    ybtd = ncread(['N',UTMlist(iutm,:),'_2010LC030_ATS_WGS.nc'],'lat');

    'fill'
    xfll_wght = get_fill_weight(xvec,xbtd);
    yfll_wght = get_fill_weight(yvec,ybtd);

    'assemble'
    for ixvec = 1:size(btd_mpd,2)
        if mod(ixvec,1000) == 0
            [ixvec size(btd_mpd,2)]
        end
        for iyvec = 1:size(btd_mpd,1)
            if xfll_wght{ixvec}(1,1) ~=0 && yfll_wght{iyvec}(1,1) ~=0
                
                BTD_sub = BTD(xfll_wght{ixvec}(:,1),yfll_wght{iyvec}(:,1));
                btd_mpd_val = yfll_wght{iyvec}(:,2).'*(BTD_sub.'*xfll_wght{ixvec}(:,2));
                btd_mpd(iyvec, ixvec) = btd_mpd_val;
                
    %             if sum(sum(BTD_sub)) > 0
    %                 imagesc(BTD_sub)
    %                 drawnow
    %             end
            end
        end
    end


%     pcolorjw(xvec,yvec,btd_mpd)
%     drawnow
%     pause
    
    UTMlist(iutm,:)
    [xbtd(1) ybtd(1)]
end
btd_mpd_0 = btd_mpd;

%%







% 
% cck = 0
% ind=0
% for izvec = 1:length(xfll_wght)
%    % [izvec xfll_wght{izvec}]
%     
%     if xfll_wght{izvec}(1,1) ~=0
%         xfll_wght{izvec}
%         cck = cck + xfll_wght{izvec}(end,2) + xfll_wght{izvec+1}(1,2) - 1
%         [izvec xfll_wght{izvec}(1,2) xfll_wght{izvec}(end,2) cck]
%         ind = izvec
%     end
% end







