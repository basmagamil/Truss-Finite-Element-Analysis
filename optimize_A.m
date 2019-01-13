function [About,ratiobout,stressbout,Abint,ratiobint,stressbint,Abvert,ratiobvert,stressbvert, stressb] = optimize_A(h,Arange)

for c=1:length(Arange)
Ab=Arange(c);
[ratiof, stressf, ratiob, stressb]=analyze_structure( h,Ab );
%elements 6 7 14 15 16 19
if (ratiob(6-5,1) >= 0.999) && (ratiob(6-5,1) <= 1)
    About=Ab;
    ratiobout=ratiob;
    stressbout=stressb(6-5,1);
end
%elements 8 9 12 13
if (ratiob(8-5,1) >= 0.99) && (ratiob(8-5,1) <= 1)
    Abint=Ab;
    ratiobint=ratiob;
    stressbint=stressb(8-5,1);
end
%elements 17 18
if (ratiob(17-5,1) >= 0.999) && (ratiob(17-5,1) <= 1)
    Abvert=Ab;
    ratiobvert=ratiob;
    stressbvert=stressb(17-5,1);
end

end
end

