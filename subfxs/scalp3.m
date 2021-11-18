function scp=scalp3(mri,cfg)

% Contributors: Gavin Paterson, Jan-Mathijs Schoffelen, Almudena Capilla & Joachim Gross
% University of Glasgow

tic
dim=size(mri);
se=strel([0,1,0;1,1,1;0,1,0]);
% se=strel([1,1,1;1,1,1;1,1,1]);

mri=double(mri);
mri=mri./max(max(max(mri)));
scp=mri>cfg.thr.low;
t=mri>cfg.thr.high;
scp=scp-t;

scp2=permute(scp,[1,3,2]);
scp3=permute(scp,[3,2,1]);
%erode
for i=1:cfg.ne
    scp=imerode(scp,se);
    scp2=imerode(scp2,se);
    scp3=imerode(scp3,se);
    scp=scp+permute(scp2,[1,3,2])+permute(scp3,[3,2,1]);
    scp=(scp>2);
    scp2=permute(scp,[1,3,2]);
    scp3=permute(scp,[3,2,1]);
    
end
%dilate
for i=1:cfg.ng
    scp=imdilate(scp,se);
    scp2=imdilate(scp2,se);
    scp3=imdilate(scp3,se);
    scp=scp+permute(scp2,[1,3,2])+permute(scp3,[3,2,1]);
    scp=(scp>2);
    scp2=permute(scp,[1,3,2]);
    scp3=permute(scp,[3,2,1]);
end

mri=zeros(dim);
for image=1:dim(3)
    for row=1:dim(1)
        
        ind=find(scp(row,:,image));
        maxind=max(ind);
        minind=min(ind);
        mri(row,minind,image)=1;
        mri(row,maxind,image)=1;
    end
end
for image=1:dim(3)
    for column=1:dim(2)
        ind=find(scp(:,column,image));
        maxind=max(ind);
        minind=min(ind);
        mri(minind,column,image)=1;
        mri(maxind,column,image)=1;
    end
end
for image=1:dim(2)
    for row=1:dim(1)
        ind=find(scp(row,image,:));
        maxind=max(ind);
        minind=min(ind);
        mri(row,image,minind)=1;
        mri(row,image,maxind)=1;
    end
end
scp  = find(mri);
[scp(:,1), scp(:,2), scp(:,3)] = ind2sub(dim, scp);
% size(scp)
scp=unique(scp,'rows');
toc

