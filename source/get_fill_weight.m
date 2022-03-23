function [fact] = get_fill_weight(zvec,zbtd)
%%
% zvec = yvec;
% zbtd = ybtd;


zbtd_d = zbtd(2) - zbtd(1);
zvec_d = zvec(2) - zvec(1);

zbtd_sub = cell(length(zvec),1);
for izvec = 1:length(zvec)
    
    
    zvec_sta = zvec(izvec) - zvec_d/2;
    zvec_end = zvec(izvec) + zvec_d/2;
    zbtd_sta = zbtd - zbtd_d/2;
    zbtd_end = zbtd + zbtd_d/2;
    
    %zbtd_end(zbtd_subslct(end))

    zbtd_subslct = find(...
        (zbtd_sta >= zvec_sta & zbtd_end <= zvec_end) | ...
        (zbtd_end >= zvec_end & zbtd_sta <= zvec_end) | ...
        (zbtd_end >= zvec_sta & zbtd_sta <= zvec_sta) ...
        );

%     if not(isempty(zbtd_subslct))
%         pause
%     end
%     
%     
%     zbtd_end(zbtd_subslct)
%     [zvec_sta zvec(izvec) zvec_end]

    if not(isempty(zbtd_subslct))
        if zbtd_subslct(1) == 1 && zvec_sta < zbtd_sta(zbtd_subslct(1))
           % zbtd_sub{izvec}(:,1) = [zbtd_subslct; zbtd_subslct(end)+1];
            zbtd_sub{izvec}(:,1) = zbtd_subslct;
            zbtd_sub{izvec}(1:end-1,2) = ones(length(zbtd_sub{izvec}(:,1))-1,1);
            zbtd_sub{izvec}(end,2) = ( zvec_end - zbtd_sta(zbtd_sub{izvec}(end,1)) ) / zbtd_d;
        elseif zbtd_subslct(end) == length(zbtd) && zvec_end > zbtd_end(zbtd_subslct(end))
            zbtd_sub{izvec}(:,1) = zbtd_subslct; %[zbtd_subslct(1)-1; zbtd_subslct ];
            zbtd_sub{izvec}(2:end,2) = ones(length(zbtd_sub{izvec}(:,1))-1,1);
            zbtd_sub{izvec}(1,2)   = ( zbtd_end(zbtd_sub{izvec}(1,1)) - zvec_sta) / zbtd_d;
        else
            zbtd_sub{izvec}(:,1) = zbtd_subslct;%[zbtd_subslct(1)-1; zbtd_subslct; zbtd_subslct(end)+1];
          %  zbtd_sub{izvec}(:,3) = [zbtd(zbtd_subslct(1)-1); zbtd(zbtd_subslct); zbtd(zbtd_subslct(end)+1)];
            zbtd_sub{izvec}(2:end-1,2) = ones(length(zbtd_sub{izvec}(:,1))-2,1);
            
            zbtd_sub{izvec}(1,2)   = ( zbtd_end(zbtd_sub{izvec}(1,1)) - zvec_sta) / zbtd_d;
            zbtd_sub{izvec}(end,2) = ( zvec_end - zbtd_sta(zbtd_sub{izvec}(end,1)) ) / zbtd_d;
        end
        zbtd_sub{izvec}
    else
        zbtd_sub{izvec}(1,:) =  zeros(1,2);
    end
    


%     if not(isempty(zbtd_sub{izvec-1})) && zbtd_sub{izvec - 1}(end,2) == 

end

% 
% cck = 0
% for izvec = 1:length(zvec)-1
%     if zbtd_sub{izvec}(1,1) ~=0
%         zbtd_sub{izvec}
%         cck = cck + zbtd_sub{izvec}(end,2) + zbtd_sub{izvec+1}(1,2) - 1
%         [zbtd_sub{izvec}(1,2) zbtd_sub{izvec}(end,2)]
%         pause
%     end
% end


fact = zbtd_sub;





