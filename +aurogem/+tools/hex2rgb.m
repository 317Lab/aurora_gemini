function rgb = hex2rgb(hex)
rgb = sscanf(hex,'%2x%2x%2x',[1 3])/255;
end
