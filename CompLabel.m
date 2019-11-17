function [Labels,comp_num] = CompLabel(imagex,nhood,isNeg,thresh)

% Performs component labelling for binary images (imagex)

% INPUT
% imagex: Input binary image (backgorund is 0 and obejcts are 1)
% nhood: Definition of neighbourhood (4-connectivity or 8-connectivity)
% back_ground: Background of the image (assumed to be the mode of the pixel distribution), to account for negative images as well. 
% thresh: eliminates objects below a certain size (default=300)

% OUTPUT
% Labels: a 2D array of integer labels, each integer label uniquely
% identifying a given component in the image.Bg is -1
%

if(nargin<4)
    thresh=size(imagex,1)*size(imagex,2)/100;
end

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
                         
            if(Labels(ii,jj)==0)
               
                M = [Labels(ii-1,jj-1) Labels(ii-1,jj) Labels(ii-1,jj+1); % Defines the neighbourhood matrix of labels
                 Labels(ii,jj-1) Labels(ii,jj) Labels(ii,jj+1)];
           
            switch nhood
                case 4 % Four connectivity
                    M_vect = [M(1,2), M(2,1)];
                    M_vect(M_vect == -1) = [];

                    if (~isempty(M_vect)) 

                        l1 = min(Mappings(M_vect));
                        for pp = 1:length(M_vect)
                            if(Mappings(M_vect(pp)) > l1)
                                Mappings(M_vect(pp))=l1;
                            end; % Propagate equivalence throughout matrix
                        end
                        Labels(ii,jj)=l1;
                    else
                        maxlabel=maxlabel+1;
                        Mappings(maxlabel)=maxlabel;
                        Labels(ii,jj)=maxlabel;

                    end

                case 8 % Eight connectivity
                    M_vect = [M(1,1), M(1,2), M(1,3), M(2,1)];
                    M_vect(M_vect == -1) = [];

                    if (~isempty(M_vect))
                        l1 = min(Mappings(M_vect));
                        for pp = 1:length(M_vect)
                          if(Mappings(M_vect(pp)) > l1)
                                Mappings(M_vect(pp))=l1;
                          end; % Optimise this step (maybe time consuming!)
                        end
                        Labels(ii,jj)=l1;
                    else
                        maxlabel=maxlabel+1;
                        Mappings(maxlabel)=maxlabel;
                        Labels(ii,jj)=maxlabel;

                    end
            end
          
            end
        for(j=1:size(Mappings,1))
             k=j;
            while(Mappings(k)~=k)
             k=Mappings(k);;
            end
            Mappings(j)=k;
        end
    end
end


%remapping all objects whose labels have been changed
for(j=1:size(Mappings,1))
    k=j;
    while(Mappings(k)~=k)
        k=Mappings(k);;
    end
    Mappings(j)=k;
end 
comp_num=max(Mappings);
% changing the mapped and rlated pixel objects

size_objects=zeros(length(Mappings),1);
for(ii=1:LL)
    for(jj=1:MM)
        if(Labels(ii,jj)>0)
            Labels(ii,jj)=Mappings(Labels(ii,jj));
            size_objects(Labels(ii,jj))=size_objects(Labels(ii,jj))+1;
        end
       
    end
end

obj_num=1;

% Evaluate for size
for(j=1:length(size_objects))
   if(size_objects(j)>thresh)
       size_objects(j)=obj_num;
       obj_num=obj_num+1;
   else
       size_objects(j)=0;
   end 
end

comp_num=obj_num-1;
%replacing by new object labels
for(ii=1:LL)
    for(jj=1:MM)
        if(Labels(ii,jj)>0)
            Labels(ii,jj)=size_objects(Labels(ii,jj));
        end
    end
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





% function [Labels] = CompLabel(imagex,nhood,back_ground)
% LL = size(imagex,1);
% MM = size(imagex,2);
% Labels = zeros(LL,MM);
% 
% % First iteration
% for ii = 2:LL-1
%     for jj = 2:MM-1
%          
%         if (imagex(ii,jj) ~= back_ground)
%                         
%              Labels(ii,jj) = max(Labels(:)) + 1;
%              
%              M = [Labels(ii-1,jj-1) Labels(ii-1,jj) Labels(ii-1,jj+1);
%              Labels(ii,jj-1) Labels(ii,jj) Labels(ii,jj+1)];
% 
%             switch nhood
%                 case 4
%                     M_vect = [M(1,2), M(2,1), M(2,2)];
%                     M_vect(M_vect == 0) = [];
% 
%                     if (~isempty(M_vect))
%                         l1 = min(M_vect);
%                          for pp = 1:length(M_vect)
%                              Labels(Labels == M_vect(pp)) = l1;
%                          end
%                     end
% 
%                 case 8
%                     M_vect = [M(1,1), M(1,2), M(1,3), M(2,1), M(2,2)];
%                     M_vect(M_vect == 0) = [];
% 
%                     if (~isempty(M_vect))
%                         l1 = min(M_vect);
%                          for pp = 1:length(M_vect)
%                              Labels(Labels == M_vect(pp)) = l1;
%                          end
%                                          
%                     end
%             end
%         end
%     end
% end
% 
% 
% % Second iteration
% % for ii = LL-1:-1:2
% %     for jj = MM-1:-1:2
% %          M = [Labels(ii,jj-1) Labels(ii,jj) Labels(ii,jj+1);
% %              Labels(ii+1,jj-1) Labels(ii+1,jj) Labels(ii+1,jj+1)];         
% %         
% %         switch nhood
% %             
% %             case 4
% %                 M_vect = [M(1,2), M(2,1)];
% %                 M_vect(M_vect == 0) = [];
% %                 
% %                 if (~isempty(M_vect))
% %                     l1 = min(M_vect);
% %                     l2 = max(M_vect);
% %                     Labels(ii,jj) = l1;
% %                     Labels(Labels == l2) = l1;
% %                 end
% %                 
% %             case 8
% %                 M_vect = [M(1,1), M(1,2), M(1,3), M(2,1)];
% %                 M_vect(M_vect == 0) = [];
% %                 
% %                 if (~isempty(M_vect))
% %                     l1 = min(M_vect);
% %                     Labels(ii,jj) = l1;
% %                     Labels(Labels == M(1,1)) = l1;
% %                     Labels(Labels == M(1,2)) = l1;
% %                     Labels(Labels == M(1,3)) = l1;
% %                     Labels(Labels == M(2,1)) = l1;
% %                 end
% %         end
% %     end
% end
% 
% 
%         
        
                    
                    
                    
                    
                    
                    
                
