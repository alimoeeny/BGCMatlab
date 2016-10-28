function C = CompareMeans(Ca,Cb)

x = AllV.ShapeCorr(Ca, Cb);
for j = 1:min([length(Ca.next) length(Cb.next)])
    x(j+1) = AllV.ShapeCorr(Ca.next{j}, Cb.next{j});
end
C.xcorrs = x;