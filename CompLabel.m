function [Labels] = CompLabel(imagex,nhood,back_ground)

% Performs component labelling for binary images (imagex)

% INPUT
% imagex: Input binary image
% nhood: Definition of neighbourhood (4-connectivity or 8-connectivity)
% back_ground: Background of the image (assumed to be the mode of the pixel distribution), to account for negative images as well. 

% OUTPUT
% Labels: a 2D array of integer labels, each integer label uniquely identifying a given component in the image.



LL = size(imagex,1);
MM = size(imagex,2);
Labels = zeros(LL,MM);

% First iteration
for ii = 2:LL-1
    for jj = 2:MM-1
         
        if (imagex(ii,jj) ~= back_ground)
                        
             Labels(ii,jj) = max(Labels(:)) + 1; % Adds a new label (higher than all previous labels)
             
             M = [Labels(ii-1,jj-1) Labels(ii-1,jj) Labels(ii-1,jj+1); % Defines the neighbourhood matrix of labels
             Labels(ii,jj-1) Labels(ii,jj) Labels(ii,jj+1)];

            switch nhood
                case 4 % Four connectivity
                    M_vect = [M(1,2), M(2,1), M(2,2)];
                    M_vect(M_vect == 0) = [];

                    if (~isempty(M_vect)) 
                        l1 = min(M_vect);
                         for pp = 1:length(M_vect)
                             Labels(Labels == M_vect(pp)) = l1;
                         end
                    end

                case 8 % Eight connectivity
                    M_vect = [M(1,1), M(1,2), M(1,3), M(2,1), M(2,2)];
                    M_vect(M_vect == 0) = [];

                    if (~isempty(M_vect))
                        l1 = min(M_vect);

                        for pp = 1:length(M_vect)
                            Labels(Labels == M_vect(pp)) = l1; % Optimise this step (maybe time consuming!)
                        end     

                    end
            end
        end
    end
end


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
end


        
                    
                    
                    
                    
                    
                    
                
