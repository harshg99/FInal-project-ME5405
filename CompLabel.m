function [Labels] = CompLabel(imagex,nhood,isNeg)

% Performs component labelling for binary images (imagex)

% INPUT
% imagex: Input binary image (backgorund is 0 and obejcts are 1)
% nhood: Definition of neighbourhood (4-connectivity or 8-connectivity)
% back_ground: Background of the image (assumed to be the mode of the pixel distribution), to account for negative images as well. 

% OUTPUT
% Labels: a 2D array of integer labels, each integer label uniquely
% identifying a given component in the image.Bg is -1
%


LL = size(imagex,1);
MM = size(imagex,2);
if(isNeg)
    imagex=1-imagex;
end
% background represented by -1; all objects assigned 0
Labels = imagex-1;


% First iteration

 
Mappings=[];        
maxlabel=0;

for(ii=1:LL)
    for(jj=1:MM)
            if(Labels(ii,jj)~=-1)
            
            M = [Labels(ii-1,jj-1) Labels(ii-1,jj) Labels(ii-1,jj+1); % Defines the neighbourhood matrix of labels
                 Labels(ii,jj-1) Labels(ii,jj) Labels(ii,jj+1)];
           
                       
            if(M(2,2)==0)
                maxlabel=maxlabel+1;
                Mappings(maxlabel)=maxlabel;
                Labels(ii,jj)=maxlabel;
            switch nhood
                case 4 % Four connectivity
                    M_vect = [M(1,2), M(2,1),M(2,2)];
                    M_vect(M_vect == -1) = [];

                    if (~isempty(M_vect)) 
                        l1 = min(M_vect);
                        for pp = 1:length(M_vect)
                            Mappings(M_vect(pp)) = l1; % Propagate equivalence throughout matrix
                        end
                    end

                case 8 % Eight connectivity
                    M_vect = [M(1,1), M(1,2), M(1,3), M(2,1),M(2,2)];
                    M_vect(M_vect == -1) = [];

                    if (~isempty(M_vect))
                        l1 = min(M_vect);
                        for pp = 1:length(M_vect)
                           Mappings(M_vect(pp)) = l1; % Optimise this step (maybe time consuming!)
                        end
                    end
            end
           end
        end
    end
end

%remapping all objects whose labels have been changed
for(j=1:size(Mappings,1))
    k=Mappings(j);
    while(Mappings(k)~=k)
        Mappings(j)=k;
    end
end

% changing the mapped and rlated pixel objects
for(i=1:LL)
    for(j=1:MM)
        Labels(ii,jj)=Mappings(Labels(ii,jj));
    end
end
% Second iteration: additional scans
% for ii = LL-1:-1:2
%     for jj = MM-1:-1:2
%          M = [Labels(ii,jj-1) Labels(ii,jj) Labels(ii,jj+1);
%              Labels(ii+1,jj-1) Labels(ii+1,jj) Labels(ii+1,jj+1)];         
%         
%         switch nhood
%             
%             case 4
%                 M_vect = [M(1,2), M(2,1)];
%                 M_vect(M_vect == 0) = [];
%                 
%                 if (~isempty(M_vect))
%                     l1 = min(M_vect);
%                     l2 = max(M_vect);
%                     Labels(ii,jj) = l1;
%                     Labels(Labels == l2) = l1;
%                 end
%                 
%             case 8
%                 M_vect = [M(1,1), M(1,2), M(1,3), M(2,1)];
%                 M_vect(M_vect == 0) = [];
%                 
%                 if (~isempty(M_vect))
%                     l1 = min(M_vect);
%                     Labels(ii,jj) = l1;
%                     Labels(Labels == M(1,1)) = l1;
%                     Labels(Labels == M(1,2)) = l1;
%                     Labels(Labels == M(1,3)) = l1;
%                     Labels(Labels == M(2,1)) = l1;
%                 end
%         end
%     end



        
                    
                    
                    
                    
                    
                    
                
