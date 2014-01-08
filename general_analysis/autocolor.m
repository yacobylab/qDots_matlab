function out=autocolor
  color=1;
  out = @nextcolor;

  function n = nextcolor
      colors='rgbcmk';
      n=colors(mod(color-1,end)+1);
      color=color+1;
  end
end
