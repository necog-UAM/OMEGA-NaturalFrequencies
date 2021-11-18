function outfile = findfile(filename)

a = dir;
for i = 1:length(a)
    if strfind(a(i).name,filename) > 0
        outfile = a(i).name;
    end
end
