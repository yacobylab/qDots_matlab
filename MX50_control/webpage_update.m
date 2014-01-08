% updates a text file that gets posted to the internet


%get the data that we want.
hel = getHelium;
T = 

fid = fopen('z:/qdots/notes/index.html', 'w');
fprintf(fid, 'Helium = %.0f%%\n Cold spot = %.2fK\n', hel, T);
fclose(fid);

% upload to webpage

f = ftp('mx50.hoopy.org', 'yacobylabmx50@hoopy.org', 'ups12345');
mput(f, 'z:/qdots/notes/index.html');
close(f); clear f; 