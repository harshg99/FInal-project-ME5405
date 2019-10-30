function [acc]= correl(m1,m2)
   
    % correlation of matrix m1 with marix m2
    s1_1=size(m1,1);
    s1_2=size(m1,2);
    s2_1=size(m2,1);
    s2_2=size(m2,2);
    
    if(s1_1~=s2_1 || s2_1~=s2_2)
        acc=0;
        return
    end
    acc=0.0;
    for(i=1:s1_1)
        for(j=1:s1_2)
            acc=acc+m2(i,j)*m1(i,j);
        end
    end
    %acc=ceil(acc);       
end