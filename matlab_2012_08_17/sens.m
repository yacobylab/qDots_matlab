v=[4 5 6 7 8 9];
f=[117.4 124.2 131.4 139.1 147.3 155.9].*1e6;
s=[6.75 7.09 7.59 8.1 8.53 8.82].*1e6;

%vave=6.43; vin=linspace(5.5,7.5,100); 
vave=6.5; vin=6.8; 
wave=interp1(v,f,vave);
w=interp1(v,f,vin); 
save=interp1(v,s,vave); 

wdiff=(w-wave)-save.*(vin-vave)

