function [h,ratiof,stressf] = optimize_h(hrange,Ab)

for z=1:length(hrange)
h=hrange(z);
[ratiof, stressf, ratiob, stressb]=analyze_structure( h,Ab );
if (max(ratiof) >= 0.999) && (max(ratiof) <= 1)
     break
end

end

end

