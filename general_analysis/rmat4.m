function m=rmat4(axis, angle)
% function m=rmat4(axis, angle)
  axis=axis/norm(axis);
  ct=cos(angle); st=sin(angle);
  x=axis(1); y=axis(2); z=axis(3);
  m=[1,0,0,0; ...
      0,ct+x*x*(1-ct), x*y*(1-ct)-z*st, x*z*(1-ct)+y*st; ...
     0,y*x*(1-ct)+z*st, ct+y*y*(1-ct), y*z*(1-ct)-x*st; ...
     0,z*x*(1-ct)-y*st, z*y*(1-ct)+x*st, ct+z*z*(1-ct)];
end


