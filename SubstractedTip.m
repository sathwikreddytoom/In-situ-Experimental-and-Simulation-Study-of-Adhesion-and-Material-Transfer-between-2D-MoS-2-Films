%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program... 
% (1) opens upper tip file         (uppertip.asc)
% (2) opens lowere tip file        (lowertip.asc)
% (3) interpolates both tip data   (f_A_inter, f_B_inter)
% (4) substracts the tip data      (C)
% (5) saves substracted tip data   (SubstractedTip.asc)
% (6) saves the graphs             (Substracted.png)
%
% The value of h_endA means how deep you obtain the tip data
% the unit is nm
h_endA = 7;
% Open the file of upper tip
profile_filename1 = 'uppertip.asc';
% Open the file of lower tip
profile_filename2 = 'lowertip.asc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














% Read the data
A=readmatrix(profile_filename1,'FileType','text');
B=readmatrix(profile_filename2,'FileType','text');

% Test plot original tip dat
subplot(2,2,1);
plot(A(:,1),A(:,2),'Bo')
hold on;
plot(B(:,1),B(:,2),'Ro')
axis equal
grid on;
xlabel('x[nm]');
ylabel('y[nm]');
legend('UpperTip (RawData)','LowereTip (RawData)');
hold off;

% Move the top of the tip to zero
[C I] = min([A(:,2)])
A(:,1) = A(:,1) - A(I,1);
A(:,2) = A(:,2) - C;
% Get the scaen area [nm]
h_startA = 0;
h_endA = h_endA;

% Move the top of the tip to zero
[C I] = max([B(:,2)])
B(:,1) = B(:,1) - B(I,1);
B(:,2) = B(:,2) - C;
% Get the scaen area
h_startB = 0;
h_endB = h_endA;

% Test plot after the movement
subplot(2,2,2);
plot(A(:,1),A(:,2),'Bo')
hold on;
plot(B(:,1),B(:,2),'Ro')
grid on;
axis equal
xlabel('x[nm]');
ylabel('y[nm]');

% Get the scan area
Amin=1000;
for i=1:length(A(:,1))
    if(A(i,2) <= h_endA) && (A(i,1) <= Amin)
        Amin = A(i,1)
        Numb = i
    end
end
z_startA = Amin
z_startA_num = Numb

Amax=-1000;
for i=1:length(A(:,1))
    if(A(i,2) <= h_endA) && (A(i,1) >= Amax)
        Amax = A(i,1)
        Numb = i
    end
end
z_endA = Amax
z_endA_num = Numb

% Get z values to plot, zplot is the number of the plot
zplot = 500;
z_interpolation_points_nmA = linspace(z_startA, z_endA, zplot);
z_interpolation_points_nmA = transpose(z_interpolation_points_nmA);

% Interpolate A and B
f_A_inter = interp1(A(:,1),A(:,2),z_interpolation_points_nmA,'linear');
disp(length(B(:,1)));
disp(length(B(:,2)));
f_B_inter = interp1(B(:,1),B(:,2),z_interpolation_points_nmA,'linear');
plot(z_interpolation_points_nmA,f_A_inter,'B.');
plot(z_interpolation_points_nmA,f_B_inter,'R.');
legend('UpperTip (Moved)','LowereTip (Moved)','Interpolated UpperTip','Interpolated LowerTip');
hold off;

% A substract from B leaves C
subplot(2,2,[3 4]);
C=f_A_inter - f_B_inter;
plot(z_interpolation_points_nmA,C,'G.');
hold on;
plot(z_interpolation_points_nmA,0,'K.');
grid on;
axis equal
xlabel('x[nm]');
ylabel('y[nm]');
legend('Substracted Tip');
hold off;

% Save the graphs as png data
saveas(gcf,'SubstractedTip','png');

% Write the Sustracted tip data
mergedC = horzcat(z_interpolation_points_nmA,C);
writematrix(mergedC,'SubstractedTip.asc','FileType','text','Delimiter','tab');





