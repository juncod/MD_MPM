function visualization(mpData, mesh)
nmp = length(mpData);
pos = zeros(nmp,2);

for i=1:nmp
    pos(i,1,1) = mpData(i).mpC(1);
    pos(i,2,1) = mpData(i).mpC(2);
    temp(i,1) = mpData(i).temp;
end
etpl = mesh.etpl;
coord = mesh.coord;

figure
scatter(pos(:,1),pos(:,2),20,temp,'s','filled'); grid on; grid minor; axis equal;
hold on;

patch('Faces',etpl,'Vertices',coord,'FaceColor','none');
axis tight
end