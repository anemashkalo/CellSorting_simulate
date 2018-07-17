function [checked_data]= stay_within_circle(new_data, radius_pxl,old_data)

coord_tmp =[];

%make sure cells don't go out of the circle 
for jj=1:size(new_data,1)
    coord_tmp(jj) = (new_data(jj,1).*new_data(jj,1) + new_data(jj,2).*new_data(jj,2));
   % disp(power(coord_tmp(jj),0.5))
    if power(coord_tmp(jj),0.5)> radius_pxl
       %disp('outside');
       new_data(jj,1:2)=old_data(jj,1:2);
       coord_tmp(jj) = (new_data(jj,1).*new_data(jj,1) + new_data(jj,2).*new_data(jj,2));
       %disp(power(coord_tmp(jj),0.5))
    end
end
checked_data=new_data;
end