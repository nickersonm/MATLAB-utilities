%% Helper function to generate circular planar array
% centers = arrayUCA(Nd)
% Generates a circular planar array with unit spacing, Nd diameter
function centers = arrayCPA(Nd)
    centers = (1:Nd); centers = centers-mean(centers);
    centers = unique(nchoosek([centers centers], 2), 'rows');
    centers = centers(centers(:,1).^2 + centers(:,2).^2 <= (Nd/2)^2,:);
end
