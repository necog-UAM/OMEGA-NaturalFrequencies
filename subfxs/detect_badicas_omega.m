function [badicas] = detect_badicas_omega (cfg, data, comp)

% Visual inspection of ICs
% 
% Left mouse click to indicate a badica 
% Arrow keys to move forward and backward
% The list of badicas will be automatically saved in the current folder

if isfield(cfg, 'layout') 
    lay = cfg.layout;
else
    lay = 'CTF275.lay';     
end

if isfield(cfg, 'scale')
    sc = cfg.scale;  
end

cfgx = cfg;

figure
count = 1;
badicas = [];
davg = ft_timelockanalysis([],data);
ic = 1;

while ic <= size(comp.trial{1},1)       
    
    d = comp.trial{1}(ic,:);
    tm = comp.time{1};
    
    % 1 ------------- complete time series  -------------%
    
    subplot(6,2,[1 3])
    
    plot(tm,d)
    set(gca,'XLim',[0 tm(end)])
    set(gca,'YLim',[-sc sc])
    
    % 2 ------------- topography -------------%
    
    subplot(6,2,[2 4])
        
    topo_ica = comp.topo(:,ic);
    
    data2 = davg;
    data2.avg = repmat(topo_ica,[1 size(davg.avg,2)]);
    
    cfg = [];
    cfg.layout = lay;
    cfg.figure = 'gca';

    ft_topoplotER(cfg,data2)    
    title(['IC #  ' num2str(ic)],'FontWeight','bold','FontSize',16)
%     colorbar
   
    % 3 -------------  single trials  --------------%
        
    t = 1500;    
    d([1:t+1,end-t-1:end])=0;
    [d2,I] = sort(abs(d),'descend');
 
    for j=1:4
        samp = I(j);
        for i=1:500
            pos = find(I==samp-i);
            I(pos)=NaN;
            pos = find(I==samp+i);
            I(pos)=NaN;
        end
        d2=d2(~isnan(I));
        I=I(~isnan(I));
        subplot(6,2,[2*j+3 2*j+4])

        try
            plot(tm(1):tm+2*t,d(I(j)-t:I(j)+t))
        catch
            j=j+1;
            plot(tm(1):tm+2*t,d(I(j)-t:I(j)+t))
        end
        set(gca,'YLim',[-sc sc],'XLim',[0 t*2])
    end  
 
    keydown = waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    if (keydown == 0)
        badicas(count)=ic;
        count=count+1;
    elseif value == 28
         disp(value);
         ic = ic-2;
    end
    ic = ic+1; 

end

save badicas badicas

close(clf)



