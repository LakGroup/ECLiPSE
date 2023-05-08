% FractalDimension.m

% Copyright Lewandowski, Z. and Beyenal, H.,  Center for Biofilm Engineering,
% Montana State University and the School of Chemical Engineering and Bioengineering,
% Washington State University.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met: 

% 1. If the program is modified, redistributions must include a notice
% indicating that the redistributed program is not identical to the 
% software distributed by Lewandowski, Z., and Beyenal, H.
% Fundamentals of Biofilm Research, 2007.

% 2. All advertising materials mentioning features or use of this software 
% must display the following acknowledgment: This product includes software 
%developed by Lewandowski, Z., and Beyenal, H.
% Fundamentals of Biofilm Research, 2007.

% We also request that use of this software be cited in publications as 
% Lewandowski, Z., and Beyenal, H.
% Fundamentals of Biofilm Research, 2007, CRC press.  

% The software is provided "AS-IS" and without warranty of any kind, 
% express, implied or otherwise, including, without limitation, any 
% warranty of merchantability  or fitness for a particular purpose. 
% In no event shall the Center for Biofilm Engineering,
% Montana State University or the School of Chemical and Bioengineering,
% Washington State University or the authors be liable for any special, incidental, indirect 
% or consequential damages of any kind, or any damages whatsoever resulting 
% from loss of use, data or profits, whether or not advised of the 
% possibility of damage, and on any theory of liability, arising out of 
% or in connection with the use or performance of this software. 

% This code was written using MATLAB 7 (MathWorks,www.mathworks.com) and 
% may be subject to certain additional restrictions as a result. 



 
function [FD, log_perimeter, log_diameter]=MinkowskiSausageDim(XX) 
% figure;
BW=XX;
 
Dmax=17; % Maximum disk diameter  
 
Perimeter_mark=bwperim(BW,4);     % mark perimeter for dilation  
d=size(BW);   % calculate size of image 
% subplot(2,2,1), imshow(Perimeter_mark);
% title('Marked perimeter');
 
% trim image border pixels
Perimeter_mark(:,1)=[];  
d2=d(2)-1;  
Perimeter_mark(:,d2)=[]; 
 
Perimeter_mark(1,:)=[];  
d1=d(1)-1;  
Perimeter_mark(d1,:)=[]; 
 
% calculate Euclidian distance to border
ED = bwdist(Perimeter_mark,'euclidean'); % Euclidean distance map  
 
% calculate area for each diameter
S=size(ED);         % get size
lengthED=S(1)*S(2); % get 1-dimensional size
 
ED1D = reshape(ED,[lengthED,1]); % transfer it into one dimensional matrix
EDS=sort(ED1D); % sort it in ascending order means ED1d(1) has maximum value  
Dmaxindex=Dmax/2+0.5;
sumD=zeros(Dmaxindex,1); % sums of each diameter
 
scan_index=1;   % initial index number of scan  
Dindex=0;       % matrix indexes for D
 
if EDS(length(EDS)) < Dmax
    Dmax = round(EDS(length(EDS)) - 1);
end
 
for D=1:2:Dmax 
    Dindex=Dindex+1;
    radius=double(D/2);
    while (EDS(scan_index)<=radius)   % scan until larger diameter
        if (EDS(scan_index)<=radius)
            sumD(Dindex)=sumD(Dindex)+1 ;  % this line only executed if the number in the cell is smaller or equal to D
        end
        scan_index=scan_index+1;    % scan in increasing direction
    end % while
end % for
 
% calculate cumulative sums
for Dindex=2:Dmaxindex    % first one does not need previous one
    sumD(Dindex)=sumD(Dindex)+sumD(Dindex-1);
end
%%%%%%%% end calculation of areas for each diameter

% from 2 to calculated index number 
for i=1: Dmaxindex-1  % Do not count fist one because log(1) is zero and it adds error
    perimeter=sumD(i+1)/(i*2+1);    % perimeter=dilated area/diameter=index*2+1
    log_perimeter(i)=log10(perimeter);
    log_diameter(i)=log10((i*2+1));
end % for

% subplot(2,2,2),  plot(log_diameter,log_perimeter); % to see the correlation
% title('Log(diameter) vs Log(perimeter)');
% xlabel('Log(diameter)');
% ylabel('Log(perimeter)');
p = polyfit(log_diameter,log_perimeter,1); %finds the slope
FD=1-p(1);   % FD= 1- slope 
% end fractal dimension


