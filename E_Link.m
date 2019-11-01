function [Label] = E_Link(nhood,ang_thold,mag_thold,Edges)
% Edge linking of the edges based on a predetermined neighbourhood, angle and magnitude threshold
% Currently not in working condition, but determines edges and labels based on closeness of gradients and angles

% INPUT
% nhood: 4 or 8 connectivity
% ang_thold: Threshold to compare angles
% mag_thold: Threshold to compare magnitudes
% Edges: Edge array (2D)

% OUTPUT
% Label: 2D array of edge labels.


% Sobel operator

G_x = [-1 -2 -1;
    0 0 0;
    1 2 1];
 
G_y = [-1 0 1;
    -2 0 2;
    -1 0 1];

kern = size(G_x,1);

% Populate the derivative matrices (\nabla f)
Xderive = zeros(size(Edges));
Yderive = zeros(size(Edges));

for ii = 2:size(Edges,1)-1
    for jj = 2:size(Edges,2)-1
        sumX = 0;
        sumY = 0;
        for kk = 1:kern^2
            i_star = floor((kk-1)/kern)+1; i_diff = i_star - (kern - 1);
            j_star = mod((kk-1),kern)+1; j_diff = j_star - (kern - 1);
            sumX = sumX + G_x(i_star,j_star)*Edges(ii-i_diff,jj-j_diff);
            sumY = sumY + G_y(i_star,j_star)*Edges(ii-i_diff,jj-j_diff);
        end
        Xderive(ii,jj) = sumX;
        Yderive(ii,jj) = sumY;
    end
end

Xderive = abs(Xderive);
Yderive = abs(Yderive);

Mag = sqrt(Xderive.^2 + Yderive.^2);
Mag = (Mag./max(Mag(:)))*255;
Angle = atan2(Yderive,Xderive);
Angle = Angle + 2*pi*(Angle < 0);

%%
LL = size(Edges,1);
MM = size(Edges,2);
Label = zeros(LL,MM); 

% Top to down
for ii = 2:LL-1
    for jj = 2:MM-1
        if (Edges(ii,jj) ~= 0)
            M = [Mag(ii-1,jj-1) Mag(ii-1,jj) Mag(ii-1,jj+1);
                Mag(ii,jj-1) Mag(ii,jj) Mag(ii,jj+1);
                Mag(ii+1,jj-1) Mag(ii+1,jj) Mag(ii,jj+1)];
            A = [Angle(ii-1,jj-1) Angle(ii-1,jj) Angle(ii-1,jj+1);
                Angle(ii,jj-1) Angle(ii,jj) Angle(ii,jj+1);
                Angle(ii+1,jj-1) Angle(ii+1,jj) Angle(ii,jj+1)];
            M = abs(Mag(ii,jj) - M);
            A = abs(Angle(ii,jj) - A);
            I = (A > pi); B = 2*I.*(pi - A);
            A = A + B;
            Label(ii,jj) = max(Label(:)) + 1;
            switch nhood
                case 4
                    M_vect = [M(1,2), M(2,1)];%,  M(2,3), M(3,2)];
                    A_vect = [A(1,2), A(2,1)];%,  A(2,3), A(3,2)];
                    i_count = [ii-1,ii];%,  ii+1,ii+1]; 
                    j_count = [jj,jj-1];%,  jj+1,jj];
                    Check = (M_vect < mag_thold).* (A_vect < ang_thold);
                    for kk = 1:length(Check)
                        if (Check(kk) ~= 0)
                            Vect(kk) = Label(i_count(kk),j_count(kk));
                        end
                    end
                    D_check = Check;
                    D_check(D_check == 0) = []; 
                    if (size(D_check) ~= 0)
                        Label(ii,jj) = min(Vect(:));
                        for kk = 1:length(Check)
                            if (Check(kk) ~= 0)
                                Label(i_count(kk),j_count(kk)) = Label(ii,jj);
                            end
                        end
                    end
                    
                case 8
                    M_vect = [M(1,1),M(1,2),M(1,3),M(2,1)]; %,M(2,3),M(3,1),M(3,2),M(3,3)];
                    A_vect = [Angle(1,1),Angle(1,2),Angle(1,3),Angle(2,1)]; %,Angle(2,3),Angle(3,1),Angle(3,2),Angle(3,3)];
                    i_count = [ii-1,ii-1,ii-1,ii]; %,ii,ii+1,ii+1,ii+1]; 
                    j_count = [jj-1,jj,jj+1,jj-1]; %,jj+1,jj-1,jj,jj+1];
                    Check = (M_vect < mag_thold).*(A_vect < ang_thold);
                    for kk = 1:length(Check)
                        if (Check(kk) ~= 0)
                            Vect(kk) = Label(i_count(kk),j_count(kk));
                        end
                    end
                    D_check = Check;
                    D_check(D_check == 0) = []; 
                    if (size(D_check) ~= 0)
                        Label(ii,jj) = min(Vect(:));
                        for kk = 1:length(Check)
                           if (Check(kk) ~= 0)
                               Label(i_count(kk),j_count(kk)) = Label(ii,jj);
                           end
                        end
                    end
            end 
        end
    end
end

%% Second iteration: bottom up
for ii = LL-1:-1:2
    for jj = 2:MM-1
        if (Edges(ii,jj) ~= 0)
            switch nhood
                case 4
                L1 = Label(ii,jj-1);
                L2 = Label(ii+1,jj);
                L3 = Label(ii,jj);
                Vect = [L1, L2, L3];
                Vect(Vect == 0) = [];
                if (length(DVect) >= 1)
                    Label(ii,jj) = min(Vect(:));
                    for k = 1:length(Vect)
                        Label(Label == Vect(k)) = min(Vect(:));
                    end
                end
                
                case 8
                L1 = Label(ii,jj-1);
                L2 = Label(ii+1,jj);
                L3 = Label(ii+1,jj-1);
                L4 = Label(ii,jj);
                Vect = [L1, L2, L3,L4];
                Vect(Vect == 0) = [];
                if (length(Vect) >= 1)
                    Label(ii,jj) = min(Vect(:));
                    for k = 1:length(Vect)
                        Label(Label == Vect(k)) = min(Vect(:));
                    end
                end
                
                end     
        end
    end
end
end



%% Previous code
% % Move from top to bottom
% for ii = 2:L-1
%     for jj = 2:M-1  
%         Label(ii,jj) = max(Label(:)) + 1;
%         switch nhood
%             %%      
%             case 4
%                 count = 0; count = cast(count,'int8');
%                 mag_check1 = abs(Mag(ii,jj) - Mag(ii,jj-1));
%                 mag_check2 = abs(Mag(ii,jj) - Mag(ii-1,jj));
%                 
%                 ang_check1 = abs(Angle(ii,jj) - Angle(ii,jj-1));
%                 ang_check1 = ang_check1 + 2*(ang_check1 > pi)*(pi - ang_check1);
%                 ang_check2 = abs(Angle(ii,jj) - Angle(ii-1,jj));
%                 ang_check2 = ang_check2 + 2*(ang_check2 > pi)*(pi - ang_check2);
%                 
%                 if (mag_check1 < mag_thold && ang_check1 < ang_thold)
%                    Label(ii,jj) = Label(ii,jj-1);
%                    count = 1;
%                 end
% 
%                 if (mag_check2 < mag_thold && ang_check2 < ang_thold)
%                     
%                     if (count == 1)
%                         l1 = min(Label(ii,jj-1),Label(ii-1,jj));
%                         l2 = max(Label(ii,jj-1),Label(ii-1,jj));
%                         Label(ii,jj) = l1;
%                         Label(Label == l2) = l1;
%                     else
%                         Label(ii,jj) = Label(ii-1,jj);
%                     end
%                     count = count + 1;
%                 end
%                
%                 %%
%                 %%
%             case 8
%                 count = 0; count = cast(count,'int8'); 
%                 mag_check1 = abs(Mag(ii,jj) - Mag(ii,jj-1));
%                 mag_check2 = abs(Mag(ii,jj) - Mag(ii-1,jj));
%                 mag_check3 = abs(Mag(ii,jj) - Mag(ii-1,jj-1));
% 
%                 ang_check1 = abs(Angle(ii,jj) - Angle(ii,jj-1));
%                 ang_check1 = ang_check1 + 2*(ang_check1 > pi)*(pi - ang_check1);
%                 ang_check2 = abs(Angle(ii,jj) - Angle(ii-1,jj));
%                 ang_check2 = ang_check2 + 2*(ang_check2 > pi)*(pi - ang_check2);
%                 ang_check3 = abs(Angle(ii,jj) - Angle(ii-1,jj-1));
%                 ang_check3 = ang_check3 + 2*(ang_check3 > pi)*(pi - ang_check3);
% 
%                 if (mag_check1 < mag_thold && ang_check1 < ang_thold)
%                    Label(ii,jj) = Label(ii,jj-1);
%                    count = 1;
%                 end
% 
%                 if (mag_check2 < mag_thold && ang_check2 < ang_thold)
% 
%                     if (count == 1)
%                         l1 = min([Label(ii,jj-1),Label(ii-1,jj)]);
%                         l2 = max([Label(ii,jj-1),Label(ii-1,jj)]);
%                         Label(ii,jj) = l1;
%                         Label(Label == l2) = l1;
%                     else
%                         Label(ii,jj) = Label(ii-1,jj);
%                     end
%                     count = count + 1;
%                 end
%                 
%                 if (mag_check3 < mag_thold && ang_check3 < ang_thold)
%                     
%                     if (count ~= 0)
%                         l1 = min([Label(ii,jj-1),Label(ii-1,jj-1)]);
%                         l2 = max([Label(ii,jj-1),Label(ii-1,jj-1)]);
%                         Label(ii,jj) = l1;
%                         Label(Label == l2) = l1;
%                     else
%                         Label(ii,jj) = Label(ii-1,jj-1);
%                     end
%                 end
%         end
%     end
%     fprintf("Completed %i out of %i of iteration 1\n",ii,L-1); 
% end
% 
% Label1 = Label;
% figure()
% imagesc(Label);
% title(["Label"]);

% % Iteration 2 from bottom to top
% for ii = L-1:-1:2
%     for jj = 2:M-1
%         switch nhood
%             case 4
%                 L1 = Label(ii,jj-1);
%                 L2 = Label(ii+1,jj);
%                 L3 = Label(ii,jj);
%                 Vect = [L1, L2, L3];
%                 Vect(Vect == 0) = [];
%                 if (length(Vect) >= 1)
%                     Label(ii,jj) = min(Vect(:));
%                     for k = 1:length(Vect)
%                         Label(Label == Vect(k)) = min(Vect(:));
%                     end
%                 end
%                     
%            case 8
%                 L1 = Label(ii,jj-1);
%                 L2 = Label(ii+1,jj);
%                 L3 = Label(ii+1,jj-1);
%                 L4 = Label(ii,jj);
%                 Vect = [L1, L2, L3,L4];
%                 Vect(Vect == 0) = [];
%                 if (length(Vect) >= 1)
%                     Label(ii,jj) = min(Vect(:));
%                     for k = 1:length(Vect)
%                         Label(Label == Vect(k)) = min(Vect(:));
%                     end
%                 end
%         end
%     end
%     fprintf("Completed %i out of %i iteration 2\n",ii,L-1);
% end
% 
% figure()
% imagesc(Label);
% title(["Label"]);



% Remove 'noisy' edges
% noise_tol = 10;
% freq = his(Label,max(Label(:)))*size(Label,1)*size(Label,2);
% 
% for ii = 1:length(freq)
%     if (freq(ii) < noise_tol)
%         Label(Label == ii) = min(Label(:));
%     end
% end


  

                
        
