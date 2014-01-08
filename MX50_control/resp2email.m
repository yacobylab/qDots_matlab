function out = resp2email()

out = [];
cmd = 'C:\python27';
cmdstr = 'python MX50readmail.py';
cmd = ['Set PATH=',cmd, ';%PATH%&', cmdstr];
[a, s]=dos(cmd);

doit = strfind(s,'<title>Status?</title>'); 
if ~isempty(doit)
    doit = doit(1);
else
    return
end
%fprintf(s);
i_i=strfind(s,'<email>');
if ~isempty(i_i)
    i_i = i_i(find(i_i>doit,1));
else
    return
end
i_e= strfind(s(i_i:end),'</email>');
addr = s(i_i+(7:i_e-2));
mod_i = strfind(s(1:i_i),'<modified>');
if ~isempty(mod_i)
    mod_i = mod_i(end);
else
   return 
end
mod_e = strfind(s(mod_i:end),'</modified>');
mod_e = mod_e(1);
mday = s(mod_i+(10:19));
mtime = s(mod_i+(21:mod_e-3));
dtime = (now - (datenum([mday, ' ',mtime])-(4/24)))*86400; % in seconds
if dtime < 0
    fprintf('time travel not so good \n');
    return
else
    %fprintf('reply to msg sent %.2f seconds ago to %s \n',dtime, addr); %gmail clock 4 hrs off from EST
end
out.addr = addr;
out.dtime = dtime;

end