function tfbin = findbin(v,tf)

% find time or frequency bin within a time/freq vector (v)

vdif = v - tf;
tfbin = find(abs(vdif)==min(abs(vdif(:))));