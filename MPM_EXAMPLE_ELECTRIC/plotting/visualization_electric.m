function visualization_electric(mpData, mesh)
%Visualization of electric field MPM results
%--------------------------------------------------------------------------
% Plots electric potential distribution and current density vectors
% at material points overlaid on the background mesh.
%--------------------------------------------------------------------------
nmp = length(mpData);
pos = zeros(nmp,2);
phi_vals = zeros(nmp,1);
Jx = zeros(nmp,1);
Jy = zeros(nmp,1);

for i=1:nmp
    pos(i,1) = mpData(i).mpC(1);
    pos(i,2) = mpData(i).mpC(2);
    phi_vals(i) = mpData(i).phi;
    Jx(i) = mpData(i).J(1);
    Jy(i) = mpData(i).J(2);
end

etpl  = mesh.etpl;
coord = mesh.coord;

figure

% Left panel: electric potential (scalar field)
subplot(1,2,1);
scatter(pos(:,1), pos(:,2), 20, phi_vals, 's', 'filled');
grid on; grid minor; axis equal; axis tight;
colorbar; title('\phi  (Electric Potential)');
hold on;
patch('Faces',etpl,'Vertices',coord,'FaceColor','none');

% Right panel: current density magnitude with vector arrows
subplot(1,2,2);
J_mag = sqrt(Jx.^2 + Jy.^2);
scatter(pos(:,1), pos(:,2), 20, J_mag, 's', 'filled');
grid on; grid minor; axis equal; axis tight;
colorbar; title('|J|  (Current Density Magnitude)');
hold on;
step = max(1, floor(nmp/300));
quiver(pos(1:step:end,1), pos(1:step:end,2), ...
       Jx(1:step:end), Jy(1:step:end), 0.5, 'k');
patch('Faces',etpl,'Vertices',coord,'FaceColor','none');
end
