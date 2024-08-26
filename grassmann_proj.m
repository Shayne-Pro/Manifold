function dist_p = grassmann_proj(SY1,SY2)

MIN_THRESH = 1e-6;

same_flag = false;
if (nargin < 2)
    SY2 = SY1;
    same_flag = true;
end
p = size(SY1,2);


[~,~,number_sets1] = size(SY1);
[~,~,number_sets2] = size(SY2);

dist_p = zeros(number_sets2,number_sets1);

if (same_flag)
    %SY1 = SY2
    for tmpC1 = 1:number_sets1
        Y1 = SY1(:,:,tmpC1);
        for tmpC2 = tmpC1:number_sets2
            tmpMatrix = Y1'* SY2(:,:,tmpC2);
            tmpProjection_Kernel_Val = sum(sum(tmpMatrix.^2));
            if (tmpProjection_Kernel_Val < MIN_THRESH)
                tmpProjection_Kernel_Val = 0;
            elseif (tmpProjection_Kernel_Val > p)
                tmpProjection_Kernel_Val = p;
            end
            dist_p(tmpC2,tmpC1) = tmpProjection_Kernel_Val;
            dist_p(tmpC1,tmpC2) = dist_p(tmpC2,tmpC1);
        end
    end
else
    for tmpC1 = 1:number_sets1
        Y1 = SY1(:,:,tmpC1);
        for tmpC2 = 1:number_sets2
            tmpMatrix = Y1'* SY2(:,:,tmpC2);
            tmpProjection_Kernel_Val = sum(sum(tmpMatrix.^2));
            if (tmpProjection_Kernel_Val < MIN_THRESH)
                tmpProjection_Kernel_Val = 0;
            elseif (tmpProjection_Kernel_Val > p)
                tmpProjection_Kernel_Val = p;
            end
            dist_p(tmpC2,tmpC1) = tmpProjection_Kernel_Val;
        end
    end
end


end