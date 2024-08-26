% Author:
% - Mehrtash Harandi (mehrtash.harandi at gmail dot com)
%
% This file is provided without any warranty of
% fitness for any purpose. You can redistribute
% this file and/or modify it under the terms of
% the GNU General Public License (GPL) as published
% by the Free Software Foundation, either version 3
% of the License or (at your option) any later version.

% outLabel=Classify_ISCRC(qVx,M*qA,Train_label,bx,ax);

function outLabel=Classify_ISCRC(A,B,Labels,x,y)
Number_Of_Classes = max(Labels);
%
res_y = zeros(Number_Of_Classes,size(x,2));
for tmpC1 = 1:Number_Of_Classes
    classIndex = (Labels==tmpC1);
    delta_x = zeros(size(x));
    delta_x(classIndex,:) = x(classIndex,:);
    res_y(tmpC1,:) = sum((A*y - B*delta_x).^2);
end
[~,MinIndex] = min(res_y);
outLabel = MinIndex;
    

end