function emailAlert(cntl)
%function emailAlert(cntl)
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.
% Send a text message if pulse tube fails or helium is low. 
% update a webpage with the fridge pressures
% cntl can be:
% init: initializes a new alerdata. this is buggy if the instrument is
% already open. if it fails try 'clean'
% clean: deletes all the gpib objects associated with the ITC
% and then try to 'init'; should ideally add a emailAlert('clean');
% start: starts the time
% stop: stops the timer;

global alertdata;

switch cntl
    case 'init'
 
        alertdata.timer = timer('TimerFcn', @textAlert, 'Period', 600, 'ExecutionMode', 'fixedRate','Name','EmailAlertTextmsg');
        alertdata.ftptimer = timer('TimerFcn', @writewebpage, 'Period', 600, 'ExecutionMode', 'fixedRate','Name','EmailAlertWebpage');
        alertdata.ftptimer.ErrorFcn = @errfunction; %try restarting the timer with a 5 min start delay
        alertdata.timer.StopFcn = @stopperfn; %print out a stopping msg
        alertdata.ftptime.StopFcn = @stopperfn;
        %alertdata.itc = gpib('ni', 0, 24); %hard coded here and in 'clean'
        alertdata.resptimer = timer('TimerFcn', @checkEmail, 'Period', 60, 'ExecutionMode', 'fixedRate','Name','ResponseTimer');
        alertdata.itc=visa('ni', 'GPIB0::24::INSTR');
        fopen(alertdata.itc);
        alertdata.itc.EOSCharCode = 'CR';
        alertdata.itc.EOIMode = 'off';
        alertdata.itc.EOSMode = 'read';
        alertdata.lastalert = 0;
        alertdata.Tmax =6.0; % upper limit on pulse tube temperature
        alertdata.Tmin = 2.4; % lwr limit;
        alertdata.G1lim = 150; 
        alertdata.MClim = .11;
        alertdata.Helim=21; %helium limit
        alertdata.P2lim = 1e-3;

    case 'start'
        %alertdata.timer.ExecutionMode = 'fixedRate';
        start(alertdata.timer);
        start(alertdata.resptimer);
        %start(alertdata.ftptimer);
    case 'stop'
        stop(alertdata.timer);
        stop(alertdata.ftptimer);
        stop(alertdata.resptimer);
        
    case 'clean'
        s = 'GPIB0-24'; % hard coded here and in the creation of gpib opject
        ind = [];
        a = instrfind;
        for i = 1:length(a);
            if strcmp(a(i).name, s)
                ind(end+1) = i;
            end
        end
        a(ind);
        
end
end

function textAlert(timerhandle,timerevent)
global alertdata;

fprintf(alertdata.itc, '%s\r', 'R1')
%temp = fscanf(itc, '%*1c%f');
ret = fscanf(alertdata.itc);
temp = str2double(ret(2:end));
helium = getHelium;
MC = getIGH('M/C');
G1 = getIGH('G1');
P2 = getIGH('P2');


% send annoying text msg if the temperature/level/pressures are bad
if ((temp > alertdata.Tmax) || (temp < alertdata.Tmin))||(helium<alertdata.Helim) ||...
        (G1 > alertdata.G1lim) || (MC >alertdata.MClim) || P2<alertdata.P2lim && ...
        (now > alertdata.lastalert+.5/24)
    pause(10); %pause and recheck. makes this robust against faulty ox inst communication
    if ((temp > alertdata.Tmax) || (temp < alertdata.Tmin))||(helium<alertdata.Helim) ||...
        (G1 > alertdata.G1lim) || (MC >alertdata.MClim)|| P2<alertdata.P2lim ...
        && (now > alertdata.lastalert+.5/24)
    
    msg= ['Cold spot is ', ret(2:end), 'K  Helium = ',...
        num2str(helium), '%', 'G1 = ', num2str(G1), 'mBar.  MC =', num2str(MC), 'K  at ',...
        datestr(now, 'mmmm dd, yyyy HH:MM:SS AM'),'. Next alert in 30 min.'];
    send_text_message_MX50('201-693-3002','Verizon', msg); %MDS        
    %send_text_message_MX50('617-945-0308','att', msg); %OD
    send_text_message_MX50('617-921-8094','Verizon', msg); %SPH
    alertdata.lastalert=now;
    %fprintf(['the temperature at the cold spot is ', ret(2:end)]);
    end
end
if ((temp < alertdata.Tmax) && (temp > alertdata.Tmin)) && (helium>alertdata.Helim) &&...
        (G1 < alertdata.G1lim) && (MC <alertdata.MClim) && P2>alertdata.P2lim && (now < alertdata.lastalert+.5/24)
    msg= ['Things are happy!  Cold spot is ', ret(2:end),'K  Helium = ',...
        num2str(helium), '%', ' at ', datestr(now, 'mmmm dd, yyyy HH:MM:SS AM'),'.'];
    send_text_message_MX50('201-693-3002','Verizon', msg);        
    send_text_message_MX50('617-945-0308','att', msg);
    alertdata.lastalert=now-1;
    %fprintf(['the temperature at the cold spot is ', ret(2:end)]);
end

end

function writewebpage(timerhandle, timerevent)
global alertdata;
fprintf(alertdata.itc, '%s\r', 'R1')
%temp = fscanf(itc, '%*1c%f');
ret = fscanf(alertdata.itc);
temp = str2double(ret(2:end));
helium = getHelium;
MC = getIGH('M/C');
G1 = getIGH('G1');

%display helium and coldspot temp on mx50.hoopy.org
fid = fopen('z:/qdots/notes/index.html', 'w');
ctime = datestr(now);
fprintf(fid, [ctime, '     Helium = %.1f%%\n        Cold spot = %.2fK\n', ...
    'M/C = %.4fK       condenser pressure = %.2fmBar'], helium, temp, MC, G1);
fclose(fid);
f = ftp('mx50.hoopy.org', 'yacobylabmx50@hoopy.org', 'ups12345');
mput(f, 'z:/qdots/notes/index.html');
close(f); clear f; 
end

function checkEmail(timerhandle, timerevent)

addr = resp2email();
if isempty(addr)
    return
end
if isstruct(addr) && isfield(addr,'dtime') && addr.dtime < 65
    addr = addr.addr;
    global alertdata;
    fprintf(alertdata.itc, '%s\r', 'R1')
    ret = fscanf(alertdata.itc);
    temp = str2double(ret(2:end));
    helium = getHelium;
    MC = getIGH('M/C');
    G1 = getIGH('G1');
    
    ctime = datestr(now);
    msg = sprintf([ctime, '     Helium = %.1f%%\n        Cold spot = %.2fK\n', ...
        'M/C = %.4fK       condenser pressure = %.2fmBar'], helium, temp, MC, G1);
    
    setpref('Internet','E_mail','yacobylabMX50@gmail.com');
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','SMTP_Username','yacobylabMX50@gmail.com');
    setpref('Internet','SMTP_Password','ups12345');
    % The following four lines are necessary only if you are using GMail as
    % your SMTP server. Delete these lines wif you are using your own SMTP
    % server.
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');
    sendmail(addr,'Status!',msg)
end
end

function errfunction(timerhandle, timerevent)
% try to restart the timer with a 5 min start delay
global alertdata;
stop(alertdata.ftptimer);
t = get(alertdata.ftptimer);
t.StartDelay = 300;
start(alertdata.ftptimer);
end

function stopperfn(timerhandle, timerevent)
global alertdata;
Tinfo = get(timerhandle); 
fprintf('Stopping Timer %s \n', Tinfo.Name);

end