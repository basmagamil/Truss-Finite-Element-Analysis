clear all
%optimizing h
hrange=linspace(0.0001,0.01,2000);
Ab=0.5e-4; %guess
[h,ratiof,stressf]=optimize_h(hrange,Ab);
%optimizing A
Arange=linspace(0.5e-4,4e-4,2000);
[About,ratiobout,stressbout,Abint,ratiobint,stressbint,Abvert,ratiobvert,stressbvert, stressb]=optimize_A(h,Arange);
%Arranging stress vector
stress=[];
stress([1,2,3,4,5],1)=stressf;
stress([6,7,14,15,16,19],1)=stressbout;
stress([8,9,12,13],1)=stressbint;
stress([17,18],1)=stressbvert;
stress(10,1)=stressb(10-5,1);
stress(11,1)=stressb(11-5,1);
%FINAL OUTPUTS
fprintf('Optimal value for the height of the frame h is %f\n',h);
fprintf('Optimal value for the Area of bar elements number (6,7,14,15,16,19) is %f\n',About);
fprintf('Optimal value for the Area of bar elements number (8,9,12,13) is %f\n',Abint);
fprintf('Optimal value for the Area of bar elements number (17,18) is %f\n',Abvert);
elemn={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19'};
output=table(stress,'rowNames',elemn);
disp(output)



