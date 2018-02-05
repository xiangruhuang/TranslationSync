% %load data/Artsquad.mat;
% loss = 0.0;
% for e = 1:222044
%     i = I(1, e);
%     j = I(2, e);
%     assert(i <= 6514);
%     assert(j <= 6514);
%     if isnan(Rgt(:, :, i))
%        continue; 
%     end
%     if isnan(Rgt(:, :, j))
%        continue; 
%     end
% %     Rgt(:, :, i)
% %     j
% %     Rgt(:, :, j)
% %     RR(:, :, e)
%     loss = loss + norm(Rgt(:, :, i)'*Rgt(:, :, j) - RR(:,:,e));
% end
% loss/222044.0

fscanf('%d %d %f %f %f\n', );